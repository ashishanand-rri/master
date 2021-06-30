#!/opt/local/bin/python
#
# v0: Flag narrow band RFI
# v1: Include flagging for vialite RFI based on average autocorrelations 11 and 22,
#		also include flagging of outliers as was present in earlier versions of AvgRecFiles
# v2: Provide an option to correct the data for Total Efficiency
# 
#
# At end
# data is stored as pol RR, fits are stored as LL and residuals to fits are stored as RL
#
# USAGE: s3RFIflag_v2.py -i <inputfile> -o 0/DG/CH
#		 s3RFIflag_v2.py -h 
#
#

import numpy as np
import aipy as a
import matplotlib.pyplot as plt
import numpy.ma as ma
import sys, getopt
import os
import time
import itertools
import math
from math import exp, expm1, sqrt, sin, cos
from read_miriad import read_file11
from read_miriad import read_file22
import pickle



start_time_beg = time.time()
start_time = time.time()

PI = np.pi

##########################################################################

def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    # Fast and numerically precise:
    variance = np.average((values-average)**2, weights=weights)
    return (average, sqrt(variance))

###########################################################################

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

##########################################################################

def median_filter_e(y, flag, tol):
	ym      = np.ma.array(y, mask=np.logical_not(flag)) 
	dev     = ym - np.mean(ym)
	abs_dev = np.abs( dev )
	mad     = np.ma.median(abs_dev)
	std     = 1.4286 * mad 
	nflag=0
	for i in range(0, len(y)):
		if flag[i]==0:
			continue # Do nothing if point is already flagged
		else:
			if dev[i] > tol*std:
				flag[i]=0 # Flag if deviation is above threshold
				nflag += 1
	print ("Number of spectral points flagged due to large errors: ",nflag)
	return flag

##########################################################################

def median_filter_ac(freq, y, flag, istart, iwindow, tol, porder):

	PORDER = porder

	ym      = np.ma.array(y, mask=np.logical_not(flag))

	flag_segment = flag[istart:istart+iwindow]
	ym_segment = ym[istart:istart+iwindow]
	freq_segment = freq[istart:istart+iwindow]

	fit = np.zeros(len(freq_segment))

	ind_good = np.where(flag_segment==1)[0]

	fit_coeff = np.polyfit( (freq_segment[ind_good]-np.mean(freq_segment[ind_good]) ), ym_segment[ind_good], PORDER)
	fit[ind_good] = np.polyval( fit_coeff, (freq_segment[ind_good]-np.mean(freq_segment[ind_good])) )

	dev     = ym_segment - fit
	abs_dev = np.abs( dev )
	mad     = np.ma.median(abs_dev)
	std     = 1.4286 * mad 

	nflag = 0
	x1=[]
	y1=[]
	for i in range(istart, istart+iwindow):
		if flag[i]==0:
			continue # Do nothing if point is already flagged
		else:
			if abs_dev[i-istart] > tol*std:
				flag[i]=0 # Flag if deviation is above threshold
				x1.append(freq_segment[i-istart])
				y1.append(ym_segment[i-istart])
				nflag += 1

#	print "Number of points exceeding threshold: ",nflag

# 	ax1=plt.subplot(1,1,1)
# 	ax1.plot((freq_segment[ind_good]),(ym_segment[ind_good]),linestyle='-',c='b')
# 	ax1.plot((freq_segment[ind_good]),(fit[ind_good]),linestyle='-',c='g')
# 	ax1.plot(x1,y1,marker='o',linestyle=' ',c='r')
# 	plt.ylabel('Temperature')
# 	plt.xlabel('Freq MHz')
# 	plt.show()

	return flag

##########################################################################

def spike_filter(freq, y, flag, tol, norder):

	PORDER = norder

	ym      = np.ma.array(y, mask=np.logical_not(flag))
	flag_segment = flag[0:len(flag)]
	ym_segment = ym[0:len(ym)]
	freq_segment = freq[0:len(freq)]
	fit = np.zeros(len(freq_segment))

	ind_good = np.where(flag_segment==1)[0]
	fit_coeff = np.polyfit( (freq_segment[ind_good]-np.mean(freq_segment[ind_good]) ), ym_segment[ind_good], PORDER)
	fit[ind_good] = np.polyval( fit_coeff, (freq_segment[ind_good]-np.mean(freq_segment[ind_good])) )

# 	ax1=plt.subplot(1,1,1)
# 	ax1.plot((freq_segment[ind_good]),(ym_segment[ind_good]),linestyle='-')

	dev     = ym_segment - fit
	abs_dev = np.abs( dev )
	mad     = np.ma.median(abs_dev)
	std     = 1.4286 * mad 
	#print "STD of the difference between polyfit and spectrum: ",std

	nflag = 0
	x2=[]
	y2=[]
	for i in range (0,len(ind_good)-2):
		j1 = ym_segment[ind_good[i]]
		j2 = ym_segment[ind_good[i+1]]
		j3 = ym_segment[ind_good[i+2]]
		dif12=np.abs(j2-j1)
		dif23=np.abs(j3-j2)
		dif13=np.abs(j1-j3)
		if dif12 > tol*std and dif23 > tol*std and dif13 < tol*std :
			flag[ind_good[i+1]] = 0 # flag if the midpoint of the trio is a spike
			nflag += 1
			x2.append(freq_segment[ind_good[i+1]])
			y2.append(ym_segment[ind_good[i+1]])

	print ("Number of single spikes detected: ",nflag)

# 	ax1.plot((x2),(y2),marker='o',linestyle=' ',c='r')

	ym      = np.ma.array(y, mask=np.logical_not(flag))
	flag_segment = flag[0:len(flag)]
	ym_segment = ym[0:len(ym)]
	freq_segment = freq[0:len(freq)]
	fit = np.zeros(len(freq_segment))

	ind_good = np.where(flag_segment==1)[0]
	fit_coeff = np.polyfit( (freq_segment[ind_good]-np.mean(freq_segment[ind_good]) ), ym_segment[ind_good], PORDER)
	fit[ind_good] = np.polyval( fit_coeff, (freq_segment[ind_good]-np.mean(freq_segment[ind_good])) )

	dev     = ym_segment - fit
	abs_dev = np.abs( dev )
	mad     = np.ma.median(abs_dev)
	std     = 1.4286 * mad 
	#print "\nIterate: STD of the difference between polyfit and spectrum: ",std

	nflag = 0
	x3=[]
	y3=[]
	for i in range (0,len(ind_good)-3):
		j1 = ym_segment[ind_good[i]]
		j2 = ym_segment[ind_good[i+1]]
		j3 = ym_segment[ind_good[i+2]]
		j4 = ym_segment[ind_good[i+3]]
		dif12=np.abs(j2-j1)
		dif34=np.abs(j3-j4)
		dif14=np.abs(j1-j4)
		if dif12 > tol*std and dif34 > tol*std and dif14 < tol*std :
			flag[ind_good[i+1]] = 0 # flag if the midpoints of the quad is a spike
			flag[ind_good[i+2]] = 0 # flag if the midpoints of the quad is a spike
			nflag += 2
			x3.append(freq_segment[ind_good[i+1]])
			y3.append(ym_segment[ind_good[i+1]])
			x3.append(freq_segment[ind_good[i+2]])
			y3.append(ym_segment[ind_good[i+2]])

# 	ax1.plot((x3),(y3),marker='o',linestyle=' ',c='b')
# 	plt.ylabel('Temperature')
# 	plt.xlabel('Freq MHz')
# 	plt.show()

	print ("Number of spike pairs detected: ",nflag/2)

	return flag

##########################################################################

def median_filter_rbe(y, flag, tol):
	ym      = np.ma.array(y, mask=np.logical_not(flag)) 
	dev     = ym - np.mean(ym)
	abs_dev = np.abs( dev )
	mad     = np.ma.median(abs_dev)
	std     = 1.4286 * mad 
	nflag = 0
	for i in range(0, len(y)):
		if flag[i]==0:
			continue # Do nothing if point is already flagged
		else:
			if abs_dev[i] > tol*std:
				flag[i]=0 # Flag if deviation is above threshold
				nflag += 1
	print ("Number of spectral points flagged due to large data/error: ",nflag)
	return flag

##########################################################################

	
def main(argv):
	print("====================Entering s3RFIflag.py Script================================")
	############################################
	threshold  = 500
	ndrops     = 10
	DATASIZE   = 8193
	#PATH       = '/home/pi/Desktop/PRATUSH-master/system_codes/'
	PATH 	   = '/home/pi/Desktop/PRATUSH-master/system_codes/'
	infile_CHeff = '/home/pi/Desktop/PRATUSH-master/system_codes'
	infile_DGeff = '/home/pi/Desktop/PRATUSH-master/system_codes'


	# Indicies corresponding to the channel range that is processed.
	# Index 0 is channel number 1.
	# Freq = (index/8192)*250.0
	ilow	   = 1802   # 55 MHz
	ihigh 	   = 2788   # 85 MHz

	ORDER_FLAG    = 15
	ORDER_DATA    = 15
	TOL_ERR       = 3.0 
	TOL_DAT       = 3.0 

	############################################

	channel    = np.asarray(range(1, DATASIZE+1))  # 1, 2, 3....8193
	freq       = (channel-1) * (250.0/(DATASIZE-1))  # MHz

	# Converst between channel numbers and frequencies
	f_c   = lambda f : np.ceil((f * (DATASIZE-1)/250.0) + 1)
	c_f   = lambda c : 	(((c-1) * 250.0/(DATASIZE-1)))

	chisq = lambda p, x, y, e : np.sum((1.0/e**2)*(func(p,x)-y)**2)/np.sum(1.0/e**2)
	func_poly = lambda p, x    : np.polyval(p,x)

	TotEff_option = 0  # default no correction applied

	############################################

	try:
		opts, args = getopt.getopt(argv,"hi:o:",["file_in=","TotEff_option="])
	except getopt.GetoptError:
		print ('s3RFIflag_v2.py -i <inputfile> ')
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print ('s3RFIflag_v2.py -i <inputfile> -o 0/DG/CH')
			print (' Select option 0 for no TotEff correction, or DG/CH for applying one of the corrections')
			sys.exit()
		elif opt in ("-i", "--file_in"):
			file_in = arg
		elif opt in ("-o", "--TotEff_option"):
			TotEff_option = arg

	file_out   = file_in + '.RFIflag'

	print  ('\nInput file is : ',file_in)
	print  ('Output file is: ',file_out)
	print  ('TotEff option is: ',TotEff_option)

	if TotEff_option == '0' : TotEff_option  = 0
	if TotEff_option != 0 and TotEff_option != 'DG' and TotEff_option != 'CH' :
		print ("Error in user input of option")
		exit(0)
	if TotEff_option == 0 :	print ('\nWill not apply Total Efficiency correction to data')
	if TotEff_option == 'DG' : print ('\nWill apply DGHalli total Efficiency correction to data')
	if TotEff_option == 'CH' : print ('\nWill apply Chimale Total Efficiency correction to data')

	##########################################################################

	uvi = a.miriad.UV(PATH + file_in)
	uvi.select('polarization', -3, 0,  include=True)  # pol=-3 selects RL, which is residues
	no_records=0
	for preamble, data in uvi.all():
		no_records+=1  # To get the number of records
	uvi.rewind()

	print ("Number of records in input file: ",no_records)

	p_sn       = []
	p_source   = []
	p_time     = []
	p_lst      = []
	p_data     = np.zeros((no_records, DATASIZE), dtype=np.complex64)  # real+imag part
	p_fit      = np.zeros((no_records, DATASIZE), dtype=np.complex64)  # real+imag part
	p_res      = np.zeros((no_records, DATASIZE), dtype=np.complex64)  # real+imag part
	p_flag     = np.zeros((no_records, DATASIZE), dtype=np.int64) # integer flags

	i=0
	for preamble, data in uvi.all():
		p_sn.append(i) # Serial number for each record; starts with 0
		p_source.append(uvi['source']) # Source list
		p_time.append(uvi['time'])
		p_lst.append((12.0/PI)*uvi['lst'])
		flags     = np.logical_not(data.mask) # Miriad flags: 1 is gooddata, Mask: False is gooddata
		p_res[i]  = data.data #Data array
		p_flag[i] = flags # Miriad type Flags array
		print("p_flag is ",p_flag[i])

		print("non zero p_flag in range 1802:2788 is",np.count_nonzero(p_flag[i,1802:2788]),"for i=",i)
		print("dimension of p_flag is",p_flag.shape)
		i        += 1

	SPEC_END  = len(p_flag)/4
	print ("Number of spectra: ",SPEC_END)  	# 	The spectra are in order SPCD 9 records, then SPCF 3 records
											#	SPCD has 11, 12, 22 in RR, LL, RL.  SPCF has 12 in RR, LL, RL

	#######################################################################

	uvr = a.miriad.UV(PATH + file_in)
	uvr.select('polarization', -1, 0,  include=True)	# pol = -1 selects RR, which is data
	i=0
	for preamble, data in uvr.all():
		p_data[i]=data.data #Data array
		i += 1

	del(uvr)
	
	#######################################################################

	uvs = a.miriad.UV(PATH + file_in)
	uvs.select('polarization', -2, 0,  include=True)    # pol = -2 selects LL, which is fits

	i=0
	for preamble, data in uvs.all():
		p_fit[i]=data.data #Data array
		i += 1

	del(uvs)

	#######################################################################
	#
	#  Create the mean spectrum
	
	nchannels = ihigh-ilow+1

	ychannel=[]
	frq = []
	for i in range (ilow,ihigh+1):
		ychannel.append(i+1)
		frq.append(c_f(i+1))

	lowf = np.min(frq)
	hghf = np.max(frq)

	print (" Processing freq range: ",lowf,hghf)
	print (" 	Channel range: ",np.min(ychannel),np.max(ychannel))
	print ("		Index range: ",ilow,ihigh)

	flg = np.ones(nchannels)

	fit_seq = []
	res_seq = []
	err_seq = []
	flag_seq = []
	data_seq = []
	nseq = 0
	for count, ff in enumerate(range(3, int(4*SPEC_END), 4)):
		print("=====ff is :",ff)
		print("count is:",count)
		print("===sum is ",np.sum(p_flag[ff][ilow:ihigh+1]))
		if np.sum(p_flag[ff-3][ilow:ihigh+1])==0 : 
			continue
		fit_seq.append(p_fit[ff-3][ilow:ihigh+1].real)
		data_seq.append(p_data[ff-3][ilow:ihigh+1].real)
		res_seq.append(p_res[ff-3][ilow:ihigh+1].real)
		err_seq.append(p_data[ff][ilow:ihigh+1].imag)
		flag_seq.append(p_flag[ff-3][ilow:ihigh+1])
		nseq += 1
	print ("Accumulated ",nseq," spectra....")
   
	flag_seq = np.array(flag_seq)
	print("shape of flag_seq is",flag_seq.shape)
	
	
	fit_seq = np.array(fit_seq)
	print("shape of fit_seq is",fit_seq.shape)
	data_seq = np.array(data_seq)
	err_seq = np.array(err_seq)
	res_seq = np.array(res_seq)
	frq = np.array(frq)

	avg_res  = np.zeros(len(ychannel))
	avg_dat  = np.zeros(len(ychannel))
	avg_fit  = np.zeros(len(ychannel))
	avg_err  = np.zeros(len(ychannel))
	flg      = np.ones (len(ychannel)) #By default all channels are good

	bad_spec = []
	for j in range(0, nseq):
		if np.sum(flag_seq[j])==0:
			bad_spec.append(j)
	good_spec = np.delete(range(0, nseq), bad_spec)

	for i in range(0, len(ychannel)): 

		# With weights:	
		ind = (np.where((flag_seq[i]==1) & (fit_seq[i]!=0.0)))[0]

		if len(ind)==0:   
			avg_dat[i]  = 0	
			avg_res[i]  = 0
			avg_err[i]  = 0
			flg[i]      = 0
			avg_fit[i]  = np.sum(fit_seq[good_spec, i]) /(len(good_spec))
			continue
							
		avg_res[i]  = np.sum(1.0/err_seq[ind, i]**2 * res_seq[ind, i])/(np.sum(1.0/err_seq[ind, i]**2))
		avg_fit[i]  = np.sum(fit_seq[good_spec, i]) /(len(good_spec))
		avg_dat[i]  = np.sum(1.0/err_seq[ind, i]**2 * data_seq[ind, i])/(np.sum(1.0/err_seq[ind, i]**2))
		avg_err[i]  = 1.0/np.sqrt(np.sum(1.0/err_seq[ind, i]**2))

	err          	       = avg_err
	dat          	       = avg_fit + avg_res

	f_flag         	       = np.zeros(len(ychannel))

	##################################################################

	print ("\nFlagging the final residuals to reject outliers ... \n")

	ind_pre         = np.where(flg==1)[0]

	c_flag          = np.polyfit((np.log10(frq[ind_pre]) - \
						np.mean(np.log10(frq[ind_pre]))), dat[ind_pre], ORDER_FLAG)
	f_flag[ind_pre] = np.polyval(c_flag, (np.log10(frq[ind_pre]) - \
						np.mean(np.log10(frq[ind_pre]))))
	res             = dat - f_flag
	flg             = median_filter_rbe(res/err, flg, TOL_DAT)
	flg             = median_filter_e(err, flg, TOL_ERR)

	ind_aft			= np.where(flg==1)[0]

	print ("\nFlagging complete!\n--- %s seconds ---" % (time.time() - start_time))
	print ("Rejected ",len(ind_pre)-len(ind_aft)," spectral points of ",len(frq)," channels.")

	##################################################################
	#
	# processing cross (11) data next

	chn_low = ilow+1
	chn_hgh = ihigh+1
	ychannel11, ydata11, yflag11, yres11, yfit11, yerror11 = \
			read_file11(file_in, PATH, DATASIZE, chn_low, chn_hgh, channel)

	##################################################################

	avg_dat11  = np.zeros(len(ychannel11))
	flg11      = np.ones (len(ychannel11)) #By default all channels are good

	##################################################################

	bad_spec11 = []
	for j in range(0, len(ydata11)):
		if np.sum(yflag11[j])==0:
			bad_spec11.append(j)

	good_spec11 = np.delete(range(0, len(ydata11)), bad_spec11)

	##################################################################

	for i in range(0, len(ychannel11)): 

		# With weights:	
		ind = ( np.where( (yflag11[:,i]==1) ))[0]

		if len(ind)==0 or flg[i]==0:   	
			avg_dat11[i]  = 0
			flg11[i]      = 0
			continue
							
		avg_dat11[i]  = np.sum(ydata11[ind, i]) / len(ind)

	##################################################################
	#
	# processing cross (22) data next

	chn_low = ilow+1
	chn_hgh = ihigh+1
	ychannel22, ydata22, yflag22, yres22, yfit22, yerror22 = \
			read_file22(file_in, PATH, DATASIZE, chn_low, chn_hgh, channel)

	##################################################################

	avg_dat22  = np.zeros(len(ychannel22))
	flg22      = np.ones (len(ychannel22)) #By default all channels are good

	##################################################################

	bad_spec22 = []
	for j in range(0, len(ydata22)):
		if np.sum(yflag22[j])==0:
			bad_spec22.append(j)

	good_spec22 = np.delete(range(0, len(ydata22)), bad_spec22)

	##################################################################

	for i in range(0, len(ychannel22)): 

		# With weights:	
		ind = ( np.where( (yflag22[:,i]==1) ))[0]

		if len(ind)==0 or flg[i]==0:   	
			avg_dat22[i]  = 0
			flg22[i]      = 0
			continue
							
		avg_dat22[i]  = np.sum(ydata22[ind, i]) / len(ind)

	##################################################################
	#
	# Average the 11 and 22 AC spectra


	x0=[]
	y0=[]
	avg_dat1122 = np.zeros(len(ychannel))
	for i in range (len(frq)):
		if flg[i]==1: 
			x0.append(ychannel[i])
			avg_dat1122[i] = (avg_dat11[i]+avg_dat22[i])/2.0
			y0.append(avg_dat1122[i])

	print ("\n Identify and reject RFI spikes in the average autocorrelation spectrum")
	ind_pre         = np.where(flg==1)[0]
	tol = 4
	norder = 15
	flg = spike_filter(frq, avg_dat1122, flg, tol, norder)
	ind_aft			= np.where(flg==1)[0]
	print ("Rejected ",len(ind_pre)-len(ind_aft)," spectral points of ",len(frq)," channels.")

	print ("\n Hampel filter based on RFI in the average autocorrelation spectrum")
	ind_pre         = np.where(flg==1)[0]
	iwindow = 50
	tol = 5.0
	porder = 5
	for istart in range (0,len(ychannel)-80-iwindow+1,iwindow/2):
		flg =  median_filter_ac(frq, avg_dat1122, flg, istart, iwindow, tol, porder)
	iwindow = 100
	tol = 5.0
	porder = 5
	for istart in range (0,len(ychannel)-80-iwindow+1,iwindow/2):
		flg =  median_filter_ac(frq, avg_dat1122, flg, istart, iwindow, tol, porder)
	ind_aft			= np.where(flg==1)[0]
	print ("Rejected ",len(ind_pre)-len(ind_aft)," spectral points of ",len(frq)," channels.")

# 	xf=[]
# 	yf=[]
# 	for i in range (len(frq)):
# 		if flg[i]==1: 
# 			xf.append(ychannel[i])
# 			yf.append(avg_dat1122[i])
# 		ax1=plt.subplot(1,1,1)
# 		ax1.plot((xf),(yf),linestyle='-')
# 		plt.ylabel('Temperature')
# 		plt.xlabel('Freq MHz')
# 		plt.show()

	# clean so-far data is in avg_res, avg_fit, dat, flg, err, frq, ychannel

	##################################################################
	#
	#  Flag narrow band RFI here
	#  RFI is defined as a zone in which the residual spectrum is 
	#  continuously positive in a specified number of channels:
	crfi = 5

	xx = []
	yy = []
	ee = []
	ff = []
	for j in range (len(ychannel)):
		ff.append(flg[j])
		xx.append(frq[j])
		yy.append(avg_fit[j] + avg_res[j])
		ee.append(err[j])
	xx = np.array(xx)
	yy = np.array(yy)
	ee = np.array(ee)

	## Fit out an 8th order polynomial
	x1 = []
	y1 = []
	e1 = []
	x1index = []
	index = 0
	for i in range (len(xx)):
		if ff[i]==1 :
			x1.append(xx[i])
			y1.append(yy[i])
			e1.append(ee[i])
			x1index.append(index)
			index += 1
		else :
			x1index.append(-1)
	e1 = np.array(e1)
	z1 = np.polyfit((x1),(y1), 8, w=1.0/e1)	
	y_fit = np.array(func_poly(z1,x1))
	y_res = (y1)-(y_fit)
	avg0,rms0 = weighted_avg_and_std(y_res, 1/e1**2)
	print ("Weighted RMS residual (K): ",rms0)

	plt.figure()
	ss = file_in
	plt.plot((x1), (y_res), color = 'red', label=ss)

	# Select and flag bands of channels wider than or equal to crfi channels

	iend = len(xx)-1

	nflagrfi=0
	for istart in range (len(xx)-crfi+1):
		if ff[istart]==0 : continue
		if ff[istart]==1 and y_res[x1index[istart]]<=0.0 : continue
		flag_channel = []
		for i in range (istart,iend+1):	
			if ff[i]==0: continue			
			if ff[i]==1 and y_res[x1index[i]]<=0.0 : break
			if ff[i]==1 and y_res[x1index[i]]>0.0 : 
				flag_channel.append(i)
				if len(flag_channel) >= crfi :
					xseg=[]
					yseg=[]
					for jj in range (len(flag_channel)):
						flg[flag_channel[jj]]=0
						xseg.append(xx[flag_channel[jj]])
						yseg.append(y_res[x1index[flag_channel[jj]]])
						nflagrfi += 1
					plt.plot((xseg), (yseg), color = 'blue')
	#print "\n Number of channels flagged as rfi over >= ",crfi," channels: ",nflagrfi
			
	#recompute the rms residual, omitting channels flagged above
	xx = []
	yy = []
	ee = []
	ff = []
	for j in range (len(ychannel)):
		ff.append(flg[j])
		xx.append(frq[j])
		yy.append(avg_fit[j] + avg_res[j])
		ee.append(err[j])
	xx = np.array(xx)
	yy = np.array(yy)
	ee = np.array(ee)

	## Fit out the same 8th order polynomial again
	x1 = []
	y1 = []
	e1 = []
	x1index = []
	index = 0
	for i in range (len(xx)):
		if ff[i]==1 :
			x1.append(xx[i])
			y1.append(yy[i])
			e1.append(ee[i])
			x1index.append(index)
			index += 1
		else :
			x1index.append(-1)
	e1 = np.array(e1)
	z1n = np.polyfit((x1),(y1), 8, w=1.0/e1)	
	y_fit = np.array(func_poly(z1,x1))
	y_fitn = np.array(func_poly(z1n,x1))
	y_res = ( (y1)-(y_fit) )
	y_resn = ( (y1)-(y_fitn) )
	y_rese = ( (y1)-(y_fit) ) / e1
	rms0 = np.std(y_rese)
	
	# flag out positive and negative excursions beyond 2 sigma.
	nflagsigma = 0
	for ichan in range (len(xx)):
		if ff[ichan] == 0 : continue
		if abs(y_rese[x1index[ichan]]) > 2.0*rms0 :  
			flg[ichan]=0
			plt.plot(xx[ichan],y_res[x1index[ichan]],'b*')
			nflagsigma += 1
	plt.xlabel('Frequency (MHz)')
	plt.ylabel('Temperature (K))')
	plt.legend(loc='best')
	plt.xlim([(lowf),(hghf)])
	plt.grid()
	plt.show()

	print ("Number of channels flagged as high standard deviation: ",nflagsigma)

	#########################################################

	# Distrbiute the flags back to original array
	for ii in range (len(ychannel)):
		if flg[ii]==0 :
			for ff in range(3, int(SPEC_END)): 
				p_flag[4*ff][ychannel[ii]-1] = 0

	#########################################################
	# This section applies corrections for Total Efficiency

	# Generate the reflection efficiency

	min1=-0.1
	max1=+0.1
	arr = np.asfarray(frq)
	min_arr = 43.5
	max_arr = 88.3
	x1  = ((max1 - min1))*((arr) - float(min_arr))/(max_arr - min_arr) + min1
	p1 = [ -9.74301785e-05,   7.09783170e+00,  -2.43572006e+02,   1.92302400e+02, \
	  -7.45475638e+02,   4.74610598e+02,  -1.94697341e+02,   5.07282362e+01, \
	  -9.93167707e+00,   1.29972342e+00,  -1.28497019e-01]
	RefEff = 10.0**(np.array(func_poly(p1,x1)))+0.055

	# DGHalli measured total efficiency

# 	model = lambda p, x: 10.0**func_poly(p[:-1],x) + p[-1]
# 	p3 = [0.000000000000000000e+00, 1.485754407259609888e-01, -4.659690216121635808e-01, \
# 		  4.374814536332932713e-01, -3.458994336327604113e+00, 7.146842844833773967e+00, \
# 		  7.220487422117392917e+00, -1.735263495307209425e+01, 6.379821227574239417e-02]
# 	x3 = np.log10(frq)
# 	TotEff_DG = model(p3, x3)

	model = lambda p, x: 10.0**func_poly(p[:-1],x) + p[-1]
	file = open(infile_DGeff,'rb')
	param_dict = pickle.load(file)
	param_tot = param_dict['param']
	xf = np.log10(frq)
	xf = np.asfarray(xf)
	TotEff_DG = model(param_tot, xf)


	# Chimale measured total efficiency

	model = lambda p, x: 10.0**func_poly(p[:-1],x) + p[-1]
	file = open(infile_CHeff,'rb')
	param_dict = pickle.load(file)
	param_tot = param_dict['param']
	xf = np.log10(frq)
	xf = np.asfarray(xf)
	TotEff_CH = model(param_tot, xf)

	if TotEff_option == 'DG' or TotEff_option == 'CH' :
		for count in range(3, int(4*SPEC_END), 4):
			index = 0
			for j in range (ilow,ihigh+1):
				if TotEff_option == 'CH' : TotEffi = 1.0/TotEff_CH[index]
				if TotEff_option == 'DG' : TotEffi = 1.0/TotEff_DG[index]
				data_real = TotEffi * p_data[count][j].real
				data_imag = TotEffi * p_data[count][j].imag
				p_data[count][j] = complex(data_real,data_imag)
				data_real = TotEffi * p_fit[count][j].real
				data_imag = TotEffi * p_fit[count][j].imag
				p_fit[count][j] = complex(data_real,data_imag)
				data_real = TotEffi * p_res[count][j].real
				data_imag = TotEffi * p_res[count][j].imag
				p_res[count][j] = complex(data_real,data_imag)
				data_real = p_data[count+3][j].real
				data_imag = TotEffi * p_data[count+3][j].imag   # scaling the error stored as SPCF imaginary
				p_data[count+3][j] = complex(data_real,data_imag)

				index += 1
	
	#########################################################

	# Writing into miriad files with updated flags
	data_data      = ma.masked_array(p_data,      mask=np.logical_not(p_flag))
	data_res       = ma.masked_array(p_res,       mask=np.logical_not(p_flag))
	data_fit       = ma.masked_array(p_fit,       mask=np.logical_not(p_flag))

	uvo = a.miriad.UV(PATH + file_out, status='new')

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
	print("====================Exiting s3RFIflag.py Script================================")
if __name__ == "__main__":
	main(sys.argv[1:])
