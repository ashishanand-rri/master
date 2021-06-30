#!/opt/local/bin/python
#
# fits using Legandre polynomials to produce 2D array of residuals in time-freq
# R17: average residuals in time only using a moving window
# R18: upgrade to average in frequency as well
# for each averaging: median filter the time-freq 2D array of residuals
#
# At end
# data is stored as pol RR, fits are stored as LL and residuals to fits are stored as RL
#
# USAGE: s3pyflag_low.py -i <inputfile>
#		 s3pyflag_low.py -h 
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

start_time = time.time()



#######
# Takes in a single record.
# lists of channels, data, weights and flags. 
# For every record, it trims to the limits of channels - low and high - specified.
########
def select_freq_1d(x1, y1, w1, flag1, low, high): 
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

	x1    = x1[i_low:i_high+1]
	y1    = y1[i_low:i_high+1]
	w1    = w1[i_low:i_high+1]
	flag1 = flag1[i_low:i_high+1]

	return x1, y1, w1, flag1, i_low, i_high
 
######
# Function called by fit_leg_poly.
# Does a chi square minimization of the data wrt legandre polynomial expansion
######
def chisq_leg(p, x, y, w, flag_call):
	ind = np.where(flag_call==1)[0]
	y_c = y[ind]
	x_c = x[ind]
	w_c = w[ind]
	chisq = np.sum( (1/(w_c)**2)*(leg.legval(x_c, p) - y_c)**2 ) / np.sum(1/w_c**2)
	return chisq

######
# Does the fmin based fit to the data segment
# Uses lists of channels (x), data (y), weights (w), flags (flag)
######
def fit_leg_poly(x, y, w, flag, POLY_ORDER):

# Successive approximation
#	p00 = np.zeros(1) # Set initial guess to zero
#	for j in range(0, POLY_ORDER+1):	
#		pa = fmin(chisq_leg, p00, args=(x, y, w, flag), \
#				ftol=1.0e-5, xtol = 1.0e-3, maxiter=3e3, maxfun=3e3, disp=False) 
#		p00 = np.concatenate((pa, [0.0])) #Generate initial guess for the next iteration

# Replaced by direct computation
	if np.sum(flag)==0:
			yfit = np.zeros(len(flag))
			yres = np.zeros(len(flag))
			return yfit, yres

	ind    =  np.where(flag==1)[0]
	pa     = leg.legfit(x[ind], y[ind], POLY_ORDER, w=1.0/w[ind]**2)

	yfit = leg.legval(x, pa) #Fit is computed on all channels, whether flagged or not	
	yres   = (y - yfit)
	return yfit, yres

#######
# Takes in channel, data, flag and error
# channel is a 1D array
# other 3 have dimension (NO_AVG X channel)
#
# Averages the 2D array into 1D, with error propagation
#######
def avg_records(x, y, flag, err, inttime): 

	sum_data = np.zeros(len(y[0]))
	err_data = np.zeros(len(y[0]))
	flag_rev = np.zeros(len(y[0]))
	avg_data = np.zeros(len(y[0]))
	int_rec  = np.zeros(len(y[0]))

	for j in range(0, len(y[0])):       # Parse along frequency
	
		if np.sum(flag[:,j]) > 0:
			flag_rev[j]=1
	
		for i in range(0, len(y)):  # Parse along time
			if flag[i][j]==1:
				sum_data[j] += (1/(err[i][j]**2))*y[i][j]	
				err_data[j] += 1.0/(err[i][j]**2)	
				int_rec[j]  += inttime[i][j]
	
		if flag_rev[j]==1:	
			avg_data[j] = sum_data[j]/err_data[j]	
			err_data[j] = 1.0/np.sqrt(err_data[j])		
	
	return x, avg_data, flag_rev, err_data, int_rec

######
# Does a Hampel filtering on residual list and returns the updated flag list and threshold
# Takes in lists of residuals after Legendre fit, list of flags and value of tolerance
######
def median_filter(y, flag, tol):
	ym      = np.ma.array(y, mask=np.logical_not(flag)) 
								# For masks, False is good point and True is bad point
							    # While for miriad, 0 is bad and 1 is good
							    # Mask is used to avoid biasing of median by already flagged points
	ym      = ym - np.ma.mean(ym)
	abs_dev = np.abs(ym)
	mad     = np.ma.median(abs_dev)
	std     = 1.4286 * mad 
	for i in range(0, len(y)):
		if flag[i]==0:
			continue # Do nothing if point is already flagged
		else:
			if ym[i] > tol*std:
				flag[i]=0 # Flag if deviation is above threshold
	return flag, tol*std

######
# Does a Hampel filtering on 2D residual list and returns the updated 2D flag list and threshold
# Takes in lists of residuals after Legendre fit, list of flags and value of tolerance
######
def median_filter2d(y, error, flag, tol):
#def median_filter2d(y, flag, tol):
	ym      = np.ma.array(y, mask=np.logical_not(flag)) 
								# For masks, False is good point and True is bad point
							    # While for miriad, 0 is bad and 1 is good
							    # Mask is used to avoid biasing of median by already flagged points
	abs_dev = np.abs(ym)
	for ii in range(0, len(flag[0])): # iterate over Number of channels
		for jj in range(0, len(flag)): # iterate over Number of time records
			if flag[jj][ii]==0:
				continue # Do nothing if point is already flagged
			else:
				if abs_dev[jj][ii] > tol*error[jj][ii]:
				#if abs_dev[jj][ii] > tol*std:
					flag[jj][ii]=0 # Flag if deviation is above threshold
	return flag, tol

######
#Rescales arr list values into range min1 to max1
######
def rescale(arr, min1, max1):

	min_arr = np.amin(arr)
	max_arr = np.amax(arr)
	arr_sc  = ((max1 - min1))*((arr) - float(min_arr))/(max_arr - min_arr) + min1
	return arr_sc	

#######
# Flag spectral residuals (res1) after spectral averaging over win_len channels
# win_len is expected to be even
# Updates flags1 list with additional flags
# Outliers are identified as tol*1.4286*median absolute deviation
#######

def flag_avg_channel(ydata_avg, yflag_avg, yerror_avg, win_len, tol):
	
	data1 = np.zeros((len(ydata_avg),(len(ydata_avg[0])-win_len+1) )) 
	flag1 = np.zeros((len(ydata_avg),(len(ydata_avg[0])-win_len+1) )) # By deafult, all points are bad!
	err1  = np.zeros((len(ydata_avg),(len(ydata_avg[0])-win_len+1) ))
	
	if (win_len%2 != 0):
		print ("In flag_avg_channel: window length requested is odd! exiting ungracefully....")
		exit(0)
	
	for i in range(0, len(ydata_avg)): #Iterate over number of records
		for j in range(0, (len(ydata_avg[0])-win_len+1)): # Iterate over channels, each j is starting point of window 
			sum_temp = 0.0
			err_temp = 0.0
			for k in range(0, win_len): #For a given window start, average good points till the width of window
				if yflag_avg[i][j+k]==1:
					sum_temp        += 1.0/(yerror_avg[i][j+k])**2 * ydata_avg[i][j+k]  		
					err_temp        += 1.0/(yerror_avg[i][j+k])**2						
					flag1[i][j]      = 1	

			if flag1[i][j]==1: #Data and error propagation
				data1[i][j] = sum_temp/err_temp
				err1[i][j]  = 1.0/np.sqrt(err_temp) 

	flag1, _ = median_filter2d(data1, err1, flag1, tol) #median filter the averaged dataset
	#flag1, _ = median_filter2d(data1, flag1, tol) #median filter the averaged dataset

	#Distribute the flags back to the input array
	for i in range(0, len(data1)):
		for j in range(0, len(data1[0])):
			if flag1[i][j]==0:
				for k in range(0, win_len):
					yflag_avg[i][j + k] = 0
	return yflag_avg
	
#################################### main starts here

def main(argv):
	DATASIZE=8193
	print("****************************************Entering sypyflag_ulow.py script==============================================")
##### variables defined here #####
	tol        = 3 # 5
	POLY_ORDER = 15
	MAX_INT_LENGTH = 200  # 30 mins. max numb of 6.20 sec time records that will be aver
##### variables defined here #####

	CL = np.array([1312])   # window begins at 40 MHz
	CH = np.array([2950])	# window ends at 90 MHz

	file_in  = ''	
	file_out = ''
	try:
		opts, args = getopt.getopt(argv,"hi:o:","ifile=")
	except getopt.GetoptError:
		print ('s3pyflag_low.py -i <inputfile>')
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print ('s3pyflag_low.py -i <inputfile>')
			sys.exit()
		elif opt in ("-i", "--ifile"):
			file_in = arg
	file_out   = file_in + 'py'

	print  ('Input file is : ',file_in)
	print  ('Output file is: ',file_out)

	uvi = a.miriad.UV('/home/pi/Desktop/PRATUSH-master/system_codes/' + file_in)
	uvo = a.miriad.UV('/home/pi/Desktop/PRATUSH-master/system_codes/' + file_out, status='new')

	print ("Input File loaded! \n")

	WIN_TOT  = len(CL) # Total number of windows
	BAND = CH[0] - CL[0] + 1 # 95 MHz 
	BAND_4 = int(np.ceil(BAND/4.0))

	#################################################

	channel = np.asarray(range(1, DATASIZE+1))  # 1, 2, 3....8193
	freq = (channel-1) * (250.0/8192.0)  # MHz

	# The first read is just to find the number of records to open appropriately sized arrays
	no_records=0
	for preamble, data in uvi.all():
		no_records+=1  # To get the number of records
	uvi.rewind()

	print ("The input file has ",no_records," records ")

	# Initialize 2D arrays to store data, fit and flags
	# These are for writing to MIRIAD

	p_sn=[]
	p_source=[]
	p_lst=[]
	p_prm=[]
	p_data = np.zeros((no_records, DATASIZE), dtype=np.complex64)  # real+imag part
	p_res  = np.zeros((no_records, DATASIZE), dtype=np.complex64)  # real+imag part residuals
	p_fit  = np.zeros((no_records, DATASIZE), dtype=np.complex64) # real+imag fit
	p_flag = np.zeros((no_records, DATASIZE), dtype=np.int64) # integer flags

	i=0
	for preamble, data in uvi.all():
	
		p_sn.append(i) # Serial number for each record; starts with 0
		p_source.append(uvi['source']) # Source list
		p_lst.append(uvi['lst']) # LST list
		p_prm.append(preamble) # Preamble list

		flags=np.logical_not(data.mask) # Miriad flags: 1 is gooddata, Mask: False is gooddata
		p_data[i]=data.data #Data array

		p_flag[i]=flags # Miriad type Flags array
		i += 1

	print ("Dataset - data and flags - loaded into 2D arrays. \n")

	#############################################################################################
	# Odd numbered channels are bad and are flagged here
	# Note: channel list starts with 1, ends with 8193

	index=[]
	for i in range(0, len(channel)):
		if np.mod(channel[i],2)!=0: 
			index.append(i)

	for i in range(0, len(p_data)):
		p_flag[i][index]=0 #Set to bad data

	#############################################################################################

	#To get the starting point right and set store_start to be the starting index
	for i in p_sn:
		if (p_source[i]=='SPCD' and p_source[i+1]=='SPCD' and p_source[i+2]=='SPCD' \
				and p_source[i+3]=='SPCF'):
			store_start=i 
			break
	print ("Start flagging from rec num. ",store_start," from where we have 3 SPCD recs followed by SPCF")

	low_f  = CL[0]
	high_f = CH[len(CH)-1]

	for i in range(store_start, 4, 4):  	# To get the indices for the required frequency range 
		ydata  = p_data[store_start].real 	#Real part of cross correlation
		ylst   = p_lst[store_start]       	#LST for each record
		yflag  = p_flag[store_start]      	#Flag for each record
		yerror = p_data[store_start+3].imag #STD for each record
		yitime = p_data[store_start+3].real #Integration time for each record, not being used right now
		_, _, _, _, ilow, ihigh = select_freq_1d(channel, ydata, yerror, yflag, low_f, high_f) 
											#Extract the indices for required freq.

	print ("Will use index range: ",ilow,ihigh," to select the channels: ", CL[0], CH[len(CH)-1])

	ychannel = channel[ilow:ihigh+1]
	ychannel = rescale(ychannel, -1, 1)
	yfit     = np.zeros((no_records, len(ychannel)), dtype=np.float64)
	yres     = np.zeros((no_records, len(ychannel)), dtype=np.float64)

	lengthfromstart = len(p_data) - store_start
	residuallength = (len(p_data) - store_start)%4
	SPEC_END   = (lengthfromstart - residuallength)//4 
		# This is the number of spectra to be processed; omit autocorrs and SPCF record
	#SPEC_END = 3

	# tt_pts stores the total number of data points
	tt_pts = np.size(p_flag[range(store_start, store_start + 4*SPEC_END, 4), ilow:ihigh+1])

	# save number of flagged points prior to flagging deviations
	no_badpoints_org = np.count_nonzero(p_flag[range(store_start, store_start + 4*SPEC_END, 4), \
					ilow:ihigh+1]==0)


	'''
		The records are in following sequence:
		a. cross  (SPCD)         [i]
			i:  Real (data)
			ii: Imag (data)
		b. auto11 (SPCD)         [i+1]
			i:  Real (data)
			ii: Imag (zero)
		c. auto22 (SPCD)	 [i+2]
			i:  Real (data)
			ii: Imag (zero)
		d. info_inttime_rmsnoise  (SPCF)	 [i+3]
			i:  Real (integration time ms)
			ii: Imag (standard deviation K)
	'''

	print ("\n Spectral Fitting of all cross real spectral data")


	for chn_low, chn_hgh in zip(CL, CH) :
		print("**********************************************************")
		print ("Fitting before any flagging")
		print ("**********************************************************")
		print ("Processing range : " + str(chn_low * 250.0/8192) + " - " + str(chn_hgh * 250.0/8192) + " MHz")
		print ("**********************************************************")
	
		# Extract the indices for required frequency range chn_low to chn_hgh
		for i in range(store_start, 4, 4):  # To get the indices for the required frequency range 
			ydata  = p_data[store_start].real #Real part of cross correlation
			yflag  = p_flag[store_start]      #Flag for each record
			yerror = p_data[store_start+3].imag #STD for each record
			_, _, _, _, iilow, iihigh = select_freq_1d(channel, ydata, yerror, yflag, chn_low, chn_hgh) 
								#Extract the indices for required freq.
		yychannel = channel[iilow:iihigh+1]
		yychannel = rescale(yychannel, -1, 1)
	
		SPC_NO = 0
		for i in range(store_start, 4*SPEC_END, 4):

			if (p_source[i]!='SPCD' or p_source[i+1]!='SPCD' or p_source[i+2]!='SPCD' or p_source[i+3]!='SPCF'):
				print ("Error in record sequence!")
				break

			ydata	 = p_data[i][iilow:iihigh+1].real
			yflag	 = p_flag[i][iilow:iihigh+1]
			yerror	 = p_data[i+3][iilow:iihigh+1].imag
	

			yfit1, yres1		    	= fit_leg_poly(yychannel, ydata, yerror, yflag, POLY_ORDER)		


			#print "at 1: ",iilow,iilow+BAND
			p_res[i][iilow:(iilow+BAND)] 	          = yres1[0:BAND]
			p_fit[i][iilow:(iilow+BAND)] 	          = yfit1[0:BAND]

#			elif chn_low == CL[0]:
#				#print "at 2: ",iilow,iilow+3*BAND_4
#				p_res[i][iilow:(iilow+3*BAND_4)] 	  = yres1[0:3*BAND_4]
#				p_fit[i][iilow:(iilow+3*BAND_4)] 	  = yfit1[0:3*BAND_4]
#
#			elif chn_low == CL[len(CL)-1]:
#				#print "at 3: ",iilow+1*BAND_4,iilow+BAND
#				p_res[i][(iilow + 1*BAND_4):(iilow+BAND)] = yres1[1*BAND_4: BAND]
#				p_fit[i][(iilow + 1*BAND_4):(iilow+BAND)] = yfit1[1*BAND_4: BAND]
#
#			else:
#				#print "at 4: ",iilow+BAND_4,iilow+3*BAND_4
#				p_res[i][(iilow+BAND_4):(iilow+3*BAND_4)] = yres1[BAND_4:3*BAND_4]
#				p_fit[i][(iilow+BAND_4):(iilow+3*BAND_4)] = yfit1[BAND_4:3*BAND_4]

			#print "Fit to spectral record number " + str(SPC_NO+1) + " of " + str(SPEC_END) + " completed!" 
	
			SPC_NO  += 1


	# Set initial window length NO_AVG for time averaging
	expo   = 0
	NO_AVG = 2**expo 

	# At most time averaging is done for MAX_INT_LENGTH records (~ 1hr) 
	# This is the outermost processing loop
	# The averaged time-freq data is meadian filtered and flags are recorded for each averaging time

	while(NO_AVG <= np.amin([SPEC_END,MAX_INT_LENGTH])): 

		specn  = (SPEC_END + 1) - NO_AVG # Number of spectra after averaging

		print ("\n Time aver with Window size : " + str(NO_AVG) + " to give " + str(specn) + " spectra ")
		
		#Temporary arrays used to pass the dataset to be averaged

		ydata1  = np.zeros((NO_AVG, len(ychannel)))   
		yflag1  = np.zeros((NO_AVG, len(ychannel)))
		yerror1 = np.zeros((NO_AVG, len(ychannel)))
		yitime1 = np.zeros((NO_AVG, len(ychannel)))
		ylst1   = np.zeros(NO_AVG)

		# Arrays for storing the averaged datasets

		ydata_avg  = np.zeros((specn, len(ychannel))) 
		yflag_avg  = np.zeros((specn, len(ychannel)))
		yerror_avg = np.zeros((specn, len(ychannel)))
		yitime_avg = np.zeros((specn, len(ychannel)))

		'''
		In the following lines of this loop, we are doing the following:

		 Every chunk of NO_AVG records is extracted, averaged and put in _avg arrays.
		 The associated error in each channel is also revised 
		 In principle, we can directly do np.mean instead of storing the data in ydata1. 
		 However, since we have to make sure that flagged points are not used in mean 
		 and have to propagate the errors, we resort to this method

		When the loop ends, we have _avg arrays where we have averaged datapoints 
		as well as associated flags and errors
		'''

		j       = 0       		# Index for writing averaged spectra
		no_spec = 0  			# Counter to denote number of spectra accumulated so far
		i       = store_start 	# Index representing first record in the set being processed

		while(i < 4*SPEC_END) : # 4 is here because the array p_data has 4 records for a given time
		
			if no_spec==0:
				f_spec = i # For tracking the first spectra of moving window

			if (no_spec < NO_AVG):
				ydata1[no_spec]    = p_res[i][ilow:ihigh+1].real
				yflag1[no_spec]    = p_flag[i][ilow:ihigh+1]
				yerror1[no_spec]   = p_data[i+3][ilow:ihigh+1].imag
				yitime1[no_spec]   = p_data[i+3][ilow:ihigh+1].real
				ylst1[no_spec]     = p_lst[i]

				no_spec = no_spec + 1
				i = i + 4
	 
			if (no_spec == NO_AVG):

				if NO_AVG == 1 : # No averaging since averaging window size is unity
					ychannel, ydata_avg[j], yflag_avg[j], yerror_avg[j], yitime_avg[j] = \
					ychannel, ydata1, yflag1, yerror1, yitime1
				else:    # Average the datasets
					ychannel, ydata_avg[j], yflag_avg[j], yerror_avg[j], yitime_avg[j] = \
					avg_records(ychannel, ydata1, yflag1, yerror1, yitime1)

				no_spec=0	
				j += 1
				i = f_spec + 4 # Start the new window from the second record of the previous window
		
		# Having got the time averaged dataset, we now do median filter in 2D time and frequency domain
	
		if np.sum(yflag_avg) > 0:
			yflag_avg, _ = median_filter2d(ydata_avg, yerror_avg, yflag_avg, tol)
			#yflag_avg, _ = median_filter2d(ydata_avg, yflag_avg, tol)

		# Finally do a 2D flagging of the data after averging in frequency over a set of windows

		for win_len in [4, 8, 16]: # 16 channel average ~ 0.5MHz bandwidth
			print ("Current window length : " + str(win_len))
			yflag_avg = flag_avg_channel(ydata_avg, yflag_avg, yerror_avg, win_len, tol)	
	
		# Updating flag table 

		for ii in range(0, len(yflag_avg[0])): # iterate over Number of channels
			for jj in range(0, len(yflag_avg)): # iterate over Number of time records
				if yflag_avg[jj][ii]==0:
					for kk in range(0, NO_AVG): 
						p_flag[4*(jj+kk)+store_start][ilow + ii] = 0

		if NO_AVG == 1:  #If the data is not averaged in time, we redo the fit after first round of median filter
			for chn_low, chn_hgh in zip(CL, CH) :

				print ("**********************************************************")
				print ("Fitting after first frequency domain flagging")
				print ("**********************************************************")
				print ("Processing range : " + str(chn_low * 250.0/8192) + " - " + str(chn_hgh * 250.0/8192) + " MHz")
				print ("**********************************************************")
			
				# Extract the indices for required frequency range chn_low to chn_hgh
				for i in range(store_start, 4, 4):  # To get the indices for the required frequency range 
					ydata  = p_data[store_start].real #Real part of cross correlation
					yflag  = p_flag[store_start]      #Flag for each record
					yerror = p_data[store_start+3].imag #STD for each record
					_, _, _, _, iilow, iihigh = select_freq_1d(channel, ydata, yerror, yflag, chn_low, chn_hgh) 
										#Extract the indices for required freq.
				yychannel = channel[iilow:iihigh+1]
				yychannel = rescale(yychannel, -1, 1)
			
				SPC_NO = 0
				for i in range(store_start, 4*SPEC_END, 4):

					ydata	 = p_data[i][iilow:iihigh+1].real
					yflag	 = p_flag[i][iilow:iihigh+1]
					yerror	 = p_data[i+3][iilow:iihigh+1].imag
			
					yfit1, yres1		    	= fit_leg_poly(yychannel, ydata, yerror, yflag, POLY_ORDER)		

					p_res[i][iilow:(iilow+BAND)] 	          = yres1[0:BAND]
					p_fit[i][iilow:(iilow+BAND)] 	          = yfit1[0:BAND]

#					elif chn_low == CL[0]:
#						p_res[i][iilow:(iilow+3*BAND_4)] 	  = yres1[0:3*BAND_4]
#						p_fit[i][iilow:(iilow+3*BAND_4)] 	  = yfit1[0:3*BAND_4]
#
#					elif chn_low == CL[len(CL)-1]:
#						p_res[i][(iilow + 1*BAND_4):(iilow+BAND)] = yres1[1*BAND_4: BAND]
#						p_fit[i][(iilow + 1*BAND_4):(iilow+BAND)] = yfit1[1*BAND_4: BAND]
#
#					else:
#						p_res[i][(iilow+BAND_4):(iilow+3*BAND_4)] = yres1[BAND_4:3*BAND_4]
#						p_fit[i][(iilow+BAND_4):(iilow+3*BAND_4)] = yfit1[BAND_4:3*BAND_4]

					#print "Fit to spectral record number " + str(SPC_NO+1) + " of " + str(SPEC_END) + " completed!" 
			
					SPC_NO  += 1

			########### Find out the spectra which have abnormal mean due to RFI #############

			print ("**********************************************************")
			print ("Beginning to reject spectra with abnormally high mean...")
			print ("**********************************************************")

			mean_t = np.zeros(SPEC_END) #To store mean of segments
			flag_t = np.ones(SPEC_END) #To store flags for each spectra, initally assumed to be 1

			# Check if the spectra is already flagged; if not, register its mean of fit		
			count = 0
			for ff in range(store_start, 4*SPEC_END, 4):
				if np.sum(p_flag[ff][ilow:ihigh+1])==0: #If all channels in the band are already flagged by previous tasks
					flag_t[count] = 0
					count  += 1
					continue
				mean_t[count] = (np.mean(p_fit[ff][ilow:ihigh+1]).real)
				count += 1
	
			nbad_org = np.count_nonzero(flag_t==0)
	
			mean_t = np.asarray(mean_t)

			win_fit  = 150 #Number of records in 1 window

			# Do a median filter on sections of length win_fit

			for ff in range(0, SPEC_END, win_fit): #Iterable is starting index of a window
				if (SPEC_END-ff) < win_fit and SPEC_END > win_fit: # If the last window has less than win_fit numbers, 
										   # Take window from the bottom
					flag_t[SPEC_END-win_fit: SPEC_END], _ = \
					median_filter(mean_t[SPEC_END-win_fit: SPEC_END], flag_t[SPEC_END-win_fit: SPEC_END], tol)
					break
		
				elif (SPEC_END-ff) < win_fit and SPEC_END <= win_fit:

					print ("Warning! Number of records is less than " + str(win_fit) + ". Might affect the statistic.")
					flag_t, _ = median_filter(mean_t, flag_t, tol)
					break

				flag_t[ff:(ff + win_fit)], _ = median_filter(mean_t[ff:(ff + win_fit)], flag_t[ff:(ff + win_fit)], tol)	

			nbad_fin = np.count_nonzero(flag_t==0)

			# Distrbiute the flags back to original array
			for ff in range(0, SPEC_END): 
				if flag_t[ff]==0:
					p_flag[4*ff+store_start][ilow:ihigh+1] = 0
	
			print ("**********************************************************")
			print ("Rejected " + str(nbad_fin - nbad_org) + " Spectra on the basis of mean")
			print ("**********************************************************")

		no_badpoints_curr = np.count_nonzero(p_flag[range(store_start, store_start + 4*SPEC_END, 4), \
						ilow:ihigh+1]==0)
		add_flag = 100 * (float(no_badpoints_curr) - float(no_badpoints_org))/float(tt_pts)
		print ("     Additional flagging so far: " + str(add_flag) + " %")

		expo += 1        # Next window for averaging
		NO_AVG = 2**expo # New length of window

		if (NO_AVG > np.amin([SPEC_END,MAX_INT_LENGTH]) and \
						2**(expo-1) < np.amin([SPEC_END,MAX_INT_LENGTH])): 
			NO_AVG=np.amin([SPEC_END,MAX_INT_LENGTH])  # one last iter with averaging time set to max

	no_badpoints_fin = np.count_nonzero(p_flag[range(store_start, store_start + 4*SPEC_END, 4), \
						ilow:ihigh+1]==0)

	add_flag = 100 * (float(no_badpoints_fin) - float(no_badpoints_org))/float(tt_pts)
	tot_flag = 100 * (float(no_badpoints_fin))/float(tt_pts)

	print ("\n Additional flagging : " + str(add_flag) + " %")
	print ("Total number of flagged points : " + str(tot_flag) + " %")

	# Writing into miriad files with updated flags
	data_data = ma.masked_array(p_data, mask=np.logical_not(p_flag))
	data_res  = ma.masked_array(p_res,  mask=np.logical_not(p_flag))
	data_fit  = ma.masked_array(p_fit,  mask=np.logical_not(p_flag))

	#np.save('/Users/EoR/Desktop/merak/p_data_tim', p_data)
	#np.save('/Users/EoR/Desktop/merak/p_res_tim',  p_res)
	#np.save('/Users/EoR/Desktop/merak/p_fit_tim',  p_fit)
	#np.save('/Users/EoR/Desktop/merak/p_flag_tim', p_flag)

	i=0
	uvi.rewind()
	uvo.init_from_uv(uvi)
	for preamble, data in uvi.all():
		uvo.copyvr(uvi)
		uvo['pol']=-1  # data in RR polarization
		uvo.write(preamble, data_data[i])
		uvo['pol']=-2  # fits in LL polarizations as a diagnostic
		uvo.write(preamble, data_fit[i])
		uvo['pol']=-3  # residuals in RL polarizations as a diagnostic
		uvo.write(preamble, data_res[i])
		i += 1
	del(uvo)
	del(uvi)
	print("Processing time: --- %s seconds ---" % (time.time() - start_time))
	print("****************************************Exiting s3pyflag_ulow.py script==============================================")
if __name__ == "__main__":
	main(sys.argv[1:])


