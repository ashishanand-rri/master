#!/opt/local/bin/python
#
# At end
# data is stored as pol RR, fits are stored as LL and residuals to fits are stored as RL
#
# USAGE: s3pymedian_v8.py -i <inputfile>
#		 s3pymedian_v8.py -h 
#
# v8: modified to read in logger data and add the second channel temps to spectra
#

import numpy as np
import aipy as a
import matplotlib.pyplot as plt
import numpy.ma as ma
import numpy.polynomial.legendre as leg
import numpy.polynomial.chebyshev as cheb
from scipy.optimize import fmin
import sys, getopt
import os
import time
import itertools
import math
from scipy.interpolate import CubicSpline

start_time = time.time()


#######
# Takes in a single record.
# lists of channels, data, weights and flags. 
# For every record, it trims to the limits of channels - low and high - specified.
########
def select_freq_1d(x1, low, high): 
	prl = 0
	prh = 0
	for i in range(0, len(x1)):
		if x1[i]<=low:
			i_low=i
			prl = 1
		if x1[i]>=high:
			i_high=i
			prh = 1
			break
	if prl==1 and prh==0:
		i_high = len(x1)-1
	if prl==0 and prh==1:
		i_low = 0	

	return i_low, i_high 
	
############################################

def mean_filter_2s(sequence, flag_table, tol):

	SPEC_END  = len(sequence)
	good_ind  = np.where(flag_table==1)[0]
	sequence  = sequence - np.mean(sequence[good_ind])
	threshold = tol * np.std(sequence[good_ind])

	for jj in range(0, SPEC_END):
		if flag_table[jj]==0:
			continue
		if (sequence[jj] > threshold) or (sequence[jj] < -threshold):
			flag_table[jj] = 0
	return flag_table

############################################

def mean_filter_1s(sequence, flag_table, tol):

	SPEC_END  = len(sequence)
	good_ind  = np.where(flag_table==1)[0]
	sequence  = sequence - np.mean(sequence[good_ind])
	threshold = tol * np.std(sequence[good_ind])

	for jj in range(0, SPEC_END):
		if flag_table[jj]==0:
			continue
		if (sequence[jj] > threshold) :
			flag_table[jj] = 0
	return flag_table
	
############################################
func = lambda p, x    : np.polyval(p, x) 

def chisq_LL (p, x, y, flag, poswt, negwt):
    ind   = np.where(flag==1)[0]
    model = func(p,x)
    res   = y[ind] - model[ind]
    sum_p = 0.0
    sum_n = 0.0
    for ii in range(0, len(res)):
        if res[ii] >= 0.0:
            sum_p = sum_p + np.abs(res[ii])
        else:
             sum_n = sum_n + np.abs(res[ii])
    return (negwt*sum_n + poswt*sum_p)/float(len(x))

############################################

def main(argv):
	print("===================Entering s3pymedian_v7.py script==================")
	############################################
	tol        = 3
	POLY_ORDER = 15
	DATASIZE   = 8193
#	PATH       = '/home/pi/Desktop/PRATUSH-master/system_codes/'
	PATH 	   = '/home/pi/Desktop/PRATUSH-master/system_codes/'
	CL         = np.array([1312])   # 40 MHz
	CH         = np.array([2950])	# 90 MHz
	file_in    = ''	
	file_out   = ''
	channel    = np.asarray(range(1, DATASIZE+1))  # 1, 2, 3....8193
	freq       = (channel-1) * (250.0/DATASIZE)  # MHz

	logger_offset = (5.5/24.0) # hours offset between UTC and IST
	############################################

	try:
		opts, args = getopt.getopt(argv,"hi:l:","ifile=")
	except getopt.GetoptError:
		print ('s3pymedian_v7.py -i <inputfile> -l <loggerfile>')
		print('Error-20-oct')
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print ('s3pymedian_v8.py -i <inputfile> -l <loggerfile>')
			sys.exit()
		elif opt in ("-i", "--ifile"):
			file_in = arg
		elif opt in ("-l", "--ifile"):
			file_logger = arg
	file_out   = file_in + '.medianpyL'

	print  ('Input file is : ',file_in)
	print  ('Output file is: ',file_out)
	print  ('Logger file is: ',file_logger)

	uvo = a.miriad.UV(PATH + file_out, status='new')

	##########################################################################
	##  Read the logger file to get REF temperatures
	##
	logtime=[]
	reftemp=[]
	platetemp=[]
	#file_logger = PATH + file_logger
	#flog = open(file_logger,'r')
	#while flog.readline():
		#line = flog.readline()
		#aa = line.split(",")
		#aa0 = aa[0].split("/")
		#if aa0[0]=='2020' :
			#ctime=[]
			#reftemp.append(float(aa[2]) + 273.15)
			#platetemp.append(float(aa[1]) + 273.15)	
			#yy = int(aa0[0])
			#mn = int(aa0[1])
			#aa1 = aa0[2].split()
			#dd = int(aa1[0])
			#aa2 = aa1[1].split(":")
			#hh = int(aa2[0])
			#mm = int(aa2[1])
			#ss = float(aa2[2])
			#logtime.append([yy,mn,dd,hh,mm,ss])

	logger_JD=[]
	logger_T=[]
	logger_PT=[]
	#for i in range (len(reftemp)):
		#logger_T.append(reftemp[i])
		#logger_PT.append(platetemp[i])
		#ctime=logtime[i]
		#year = ctime[0]
		#month = ctime[1]
		#day = ctime[2]
		#hh = ctime[3]
		#mm = ctime[4]
		#ss = ctime[5]

		#if month == 1 or month == 2:
			#yearp = year - 1
			#monthp = month + 12
		#else:
			#yearp = year
			#monthp = month

		#A = math.trunc(yearp / 100.)
		#B = 2 - A + math.trunc(A / 4.)
		
		#if yearp < 0:
			#C = math.trunc((365.25 * yearp) - 0.75)
		#else:
			#C = math.trunc(365.25 * yearp)
		
		#D = math.trunc(30.6001 * (monthp + 1))
	
		#julian_date = B + C + D + day + 1720994.5

		#julian_date += ss/(3600.0*24.0) + float(mm)/(60.0*24.0) + float(hh)/24.0
		#logger_JD.append(julian_date - logger_offset)

# 	ax1=plt.subplot(1,1,1)
# 	ax1.plot((logger_JD),(logger_T),linestyle='-',c='b')
# 	ax1.plot((logger_JD),(logger_PT),linestyle='-',c='r')
# 	plt.xlabel('Julian Date')
# 	plt.ylabel('Reference temperature (K)')
# 	plt.show()

	#minloggerJD = np.min(logger_JD)
	#maxloggerJD = np.max(logger_JD)

	#print ("Logger Julian date range : ",minloggerJD,maxloggerJD)

	#LoggerT_interpolation = CubicSpline(logger_JD, logger_T)
				
	############################################

	uvi = a.miriad.UV(PATH + file_in)
	uvi.select('polarization', -3, 0,  include=True)
	no_records=0
	for preamble, data in uvi.all():
		no_records+=1  # To get the number of records
	uvi.rewind()

	print ("Number of records in input file: ",no_records)

	p_sn       = []
	p_source   = []
	p_time     = []
	p_data     = np.zeros((no_records, DATASIZE), dtype=np.complex64)  # real+imag part
	p_fit      = np.zeros((no_records, DATASIZE), dtype=np.complex64)  # real+imag part
	p_res      = np.zeros((no_records, DATASIZE), dtype=np.complex64)  # real+imag part
	p_flag     = np.zeros((no_records, DATASIZE), dtype=np.int64) # integer flags

	i=0
	for preamble, data in uvi.all():
		p_sn.append(i) # Serial number for each record; starts with 0
		p_source.append(uvi['source']) # Source list
		p_time.append(uvi['time'])
		flags     = np.logical_not(data.mask) # Miriad flags: 1 is gooddata, Mask: False is gooddata
		p_res[i]  = data.data #Data array
		p_flag[i] = flags # Miriad type Flags array
		#print("p_flag is ",p_flag[i])

		#print("non zero p_flag in range 1311:2849 is",np.count_nonzero(p_flag[i,1311:2949]),"for i=",i)
		#print("dimension of p_flag is",p_flag.shape)
	
		i        += 1
	#print("total non zero p_flag is",np.count_nonzero(p_flag))
	#print ("Data Julian date range : ",np.min(p_time),np.max(p_time))
	#if np.min(p_time) < minloggerJD or np.max(p_time) > maxloggerJD : 
		#print ("data outside times of logger recording in file: exiting ")
		#exit(0)

	SPEC_END  = len(p_flag)/4
	print ("Number of spectra: ",int(SPEC_END))

	#######################################################################

	uvr = a.miriad.UV(PATH + file_in)
	uvr.select('polarization', -1, 0,  include=True)
	i=0
	for preamble, data in uvr.all():
		p_data[i]=data.data #Data array
		
		#print("p_data is ",p_flag[i])

		#print("non zero p_data in range 1311:2849 is",np.count_nonzero(p_data[i,1311:2949]),"for i=",i)
		#print("dimension of p_flag is",p_flag.shape)
		
		i += 1

	del(uvr)
	
	#######################################################################

	uvs = a.miriad.UV(PATH + file_in)
	uvs.select('polarization', -2, 0,  include=True)

	i=0
	for preamble, data in uvs.all():
		p_fit[i]=data.data #Data array
		i += 1

	del(uvs)

	#######################################################################
	#
	#  Add reference temperature to SPCD 12 real spectra

	ilow, ihigh = select_freq_1d(channel, CL, CH)
	#print("ilow and ihigh are " , ilow,ihigh)

	xxx=[]
	yyy=[]
	for count in range(0, int(4*SPEC_END), 4):
		ctime = p_time[count]
		#res = min(enumerate(logger_JD), key=lambda x: abs(ctime - x[1]))
		#RefTemp = (logger_T[res[0]])
		#RefTemp = LoggerT_interpolation(ctime)
		RefTemp = 300.00
		xxx.append(ctime)
		yyy.append(RefTemp)
		for j in range (ilow,ihigh+1):
			data_real = p_data[count][j].real + RefTemp
			data_imag = p_data[count][j].imag
			p_data[count][j] = complex(data_real,data_imag)
			data_real = p_fit[count][j].real + RefTemp
			data_imag = p_fit[count][j].imag
			p_fit[count][j] = complex(data_real,data_imag)

	#ax1=plt.subplot(1,1,1)
	#ax1.plot((logger_JD),(logger_T),linestyle='-',c='g')
	#ax1.plot((logger_JD),(logger_PT),linestyle='-',c='r')
	#ax1.plot((xxx),(yyy),linestyle='-',c='b')
	#plt.xlabel('Julian Date')
	#plt.ylabel('Reference temperature (K)')
	#plt.show()

	# exit(0)


	#######################################################################
	
	print ("**********************************************************")
	print ("Beginning to reject spectra with abnormally high mean...")
	print ("**********************************************************")
	
	#######################################################################
	
	time_p  = np.asarray(range(0, int(SPEC_END)))
	#print("time_p is",time_p)
	flag_t  = np.ones(int(SPEC_END)) #To store flags for each spectra, initally assumed to be 1
	#print("flag_t is",flag_t)
	#######################################################################
	
	mean_RL = np.zeros(int(SPEC_END))
	std_RL  = np.zeros(int(SPEC_END)) 

	for count, ff in enumerate(range(0, int(4*SPEC_END), 4)):
		#print("=====ff is :",ff)
		#print("count is:",count)
		#print("===sum is ",np.sum(p_flag[ff][ilow:ihigh+1]))
		if np.sum(p_flag[ff][ilow:ihigh+1])==0 : 
			#print("ff is :",ff)
			#print("count is:",count)
			flag_t[count] = 0
	#print("flag_t later is",flag_t)
	#print ("Number of spectra flagged in input file: ",np.count_nonzero(flag_t==0)," of ",len(flag_t))

	for count, ff in enumerate(range(0, int(4*SPEC_END), 4)):
		if np.sum(p_flag[ff][ilow:ihigh+1])==0 or np.abs(np.mean(p_res[ff][ilow:ihigh+1].real)) > 1e6: 
			flag_t[count] = 0
			continue

		sequence       = p_res[ff][ilow:ihigh+1].real
		flag_seq       = p_flag[ff][ilow:ihigh+1].real

		ind_good       = np.where(flag_seq==1)[0]

		mean_RL[count] = np.mean(sequence[ind_good])
		std_RL[count]  = np.std(sequence[ind_good])
	
	nbad_org   = np.count_nonzero(flag_t==0)
	print ("Rejected spectra with extremely high mean residual: ",np.count_nonzero(flag_t==0)," of ",len(flag_t))

	#######################################################################
	
	flag_t     = mean_filter_2s(mean_RL, flag_t, tol)		
	flag_t     = mean_filter_2s(mean_RL, flag_t, tol)		
	print ("Rejected spectra with abs(mean residuals) exceeding ",tol," SD: ",np.count_nonzero(flag_t==0)," of ",len(flag_t))
    
	#######################################################################

	ind        = np.where(flag_t==1)[0]
	print("ind is",ind)
	print("std_RL[ind]",std_RL[ind])
	print("time_p[ind]",time_p[ind])
	c_coeff_RL = np.polyfit(time_p[ind], std_RL[ind], 5) 
	fit_RL     = np.polyval(c_coeff_RL, time_p)
	std_RL     = std_RL - fit_RL

	#######################################################################

	flag_t 	   = mean_filter_1s(std_RL, flag_t, tol)
	print ("5-order polyfit to run of mean residuals, rejected spectra exceeding ",tol," SD: ", \
													np.count_nonzero(flag_t==0)," of ",len(flag_t))

	#######################################################################

	mean_LL    = np.zeros(int(SPEC_END)) #To store mean of segments

	for count, ff in enumerate(range(0, int(4*SPEC_END), 4)):
		if np.sum(p_flag[ff][ilow:ihigh+1])==0 or np.abs(np.mean(p_fit[ff][ilow:ihigh+1].real)) > 1e6 :
			# or (np.mean(p_fit[ff][ilow:ihigh+1].real)) < 0: 
			flag_t[count] = 0
			continue
		mean_LL[count] = (np.mean(p_fit[ff][ilow:ihigh+1].real))
	print ("Rejected spectra with extremely high mean fits: ",np.count_nonzero(flag_t==0)," of ",len(flag_t))

	#########################################################

	ind      = np.where(flag_t!=0)[0]
	p00      = np.polyfit(time_p[ind], mean_LL[ind], POLY_ORDER)
	
	#########################################################

	c_LL     = fmin(chisq_LL, p00, args=(time_p, mean_LL, flag_t, 1, 4), ftol=1.0e-06,  \
														xtol = 1.0e-6,maxiter=2e5, maxfun=4e5)  
	fit_LL   = func(c_LL, time_p)
	res_LL   = mean_LL - fit_LL	

	#########################################################

	flag_t   = mean_filter_2s(res_LL, flag_t, tol)
	print ("Weighted polyfit to run of mean of fits, rejected spectra exceeding ",tol," SD: ", \
													np.count_nonzero(flag_t==0)," of ",len(flag_t))

	#########################################################

	c_LL     = fmin(chisq_LL, c_LL, args=(time_p, mean_LL, flag_t, 1, 4), ftol=1.0e-06,  \
														xtol = 1.0e-6,maxiter=2e5, maxfun=4e5)  
	fit_LL   = func(c_LL, time_p)
	res_LL   = mean_LL - fit_LL	

	#########################################################

	flag_t   = mean_filter_2s(res_LL, flag_t, tol)
	print ("Repeat weighted polyfit to run of mean of fits, rejected spectra exceeding ",tol," SD: ", \
													np.count_nonzero(flag_t==0)," of ",len(flag_t))

	#########################################################

	c_LL     = fmin(chisq_LL, c_LL, args=(time_p, mean_LL, flag_t, 1, 1), ftol=1.0e-06,  \
														xtol = 1.0e-6,maxiter=2e5, maxfun=4e5)  
	fit_LL   = func(c_LL, time_p)
	res_LL   = mean_LL - fit_LL	

	#########################################################

	flag_t   = mean_filter_2s(res_LL, flag_t, tol)
	print ("Repeat weighted polyfit to run of mean of fits, rejected spectra exceeding ",tol," SD: ", \
													np.count_nonzero(flag_t==0)," of ",len(flag_t))

	########################################################
	# Distrbiute the flags back to original array
	for ff in range(0, int(SPEC_END)): 
		if flag_t[ff]==0:
			p_flag[4*ff][ilow:ihigh+1] = 0
	
	#########################################################
	
	nbad_fin  = np.count_nonzero(flag_t==0)
	excess    = nbad_fin - nbad_org
	per       = np.round(float(excess)/int(SPEC_END)*100, decimals=2)

	print ("**********************************************************")
	print ("Rejected " + str(nbad_fin - nbad_org) + " Spectra (" + str(per) + "%) on the basis of mean rms/fits")
	print ("**********************************************************")

	#########################################################

	# Writing into miriad files with updated flags
	data_data      = ma.masked_array(p_data,      mask=np.logical_not(p_flag))
	data_res       = ma.masked_array(p_res,       mask=np.logical_not(p_flag))
	data_fit       = ma.masked_array(p_fit,       mask=np.logical_not(p_flag))

	i=0
	uvi.rewind()
	uvo.init_from_uv(uvi)
	for preamble, data in uvi.all():
		uvo.copyvr(uvi)
		uvo['pol']=-1  # data in RR polarization
		uvo.write(preamble, data_data[i])
		uvo['pol']=-2  # Legendre fits in LL polarizations as a diagnostic
		uvo.write(preamble, data_fit[i])
		uvo['pol']=-3  # residuals in RL polarizations as a diagnostic
		uvo.write(preamble, data_res[i])
		i += 1
	del(uvo)
	del(uvi)
	print("Processing time: --- %s seconds ---" % (time.time() - start_time))
	print("===================Exiting s3pymedian_v7.py script==================")
if __name__ == "__main__":
	main(sys.argv[1:])
