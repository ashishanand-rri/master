#!/opt/local/bin/python
#############
#
# AvgRecFiles_v7.py
#
# original average_multiple_nights.py from Saurabh 03Sep2017
#
# Ravi: change to saras 3 datasize
#		read in files created by s3pymedian and average all records in each file to give a single record
#		which is written into an output file.
#
# Ravi: v3->v4 
#		use average 11 & 22 spectra to detect channels with RFI from Vialite and reject all those channels.
#
# v5->v7 omit all flagging....just average and write out spectra.  Accept only one input file.
#
##############
# USAGE: AvgRecFiles_v7.py -i <inputfile1> 
#		 AvgRecFiles_v7.py -h 
##############

import numpy.ma as ma
import numpy as np
from numpy import savetxt
import aipy as a
import matplotlib.pyplot as plt
import matplotlib
matplotlib.interactive(False)
import sys, getopt
import os
import time
from scipy import interpolate
import itertools
from read_miriad import read_file
from scipy import signal
import copy
from scipy.optimize import fmin, basinhopping

##########################################################################

start_time_beg = time.time()
start_time     = time.time()

##########################################################################

def read_files(fname):

	open_file = open(fname, 'r')
	words_lst = []
	contents  = open_file.readlines()
	for i in range(len(contents)):
		words_lst.append(contents[i].strip('\n'))
	open_file.close()
	return words_lst

##########################################################################


def main(argv):

	file_in=[]	
	file_out = ''
	try:
		opts, args = getopt.getopt(argv,"hi:o:","ifile=")
	except getopt.GetoptError:
		print ('AvgRecFiles_v7.py -i <inputfile1>')
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print ('AvgRecFiles_v7.py -i <inputfile1> ')
			sys.exit()
		elif opt in ("-i", "--ifile"):
			file_in.append(arg)
			for i in range (len(args)):
				file_in.append(args[i])
	file_out   = file_in[0] + 'avg'

	print  ('\nInput file(s) is : ',file_in)
	print  ('Output file is: ',file_out)

	nfiles = len(file_in)

	if nfiles > 1 :
		print ("This version accepts only one input file!")
		exit(0)

	########################################################################

	DATASIZE = 8193

	f_c   = lambda f : np.ceil((f * (DATASIZE-1)/250.0) + 1)
	c_f   = lambda c : 	(((c-1) * 250.0/(DATASIZE-1)))

	func  = lambda p, x       : np.polyval(p,x) 
	chisq = lambda p, x, y, e : np.sum((1.0/e**2)*(func(p,x)-y)**2)/np.sum(1.0/e**2)

	##########################################################################

	channel       = np.asarray(range(1, DATASIZE+1))  # 1, 2, 3....8194

	chn_low       = 1802 # 1638    
	chn_hgh       = 2788 # 3278 

	print ("\nWill use channels in range: ",chn_low,chn_hgh," including end points")

#	PATH          = '/Users/rsubrahm/MIRSPACE/EOR_SARAS3_DATA/SARAS3_CHIMALE_March2020/'
#	PATH 		  = '/Users/rsubrahm/MIRSPACE/EOR_SARAS3_DATA/SARAS3_DGHALLI_Jan2020/'
	PATH          = '/home/pi/Desktop/PRATUSH-master/system_codes/'	

	##########################################################################


	##########################################################################
	
	
	print("file_in[0] is ",file_in[0])
	print("PATH is ",PATH)
	print("DATASIZE is ",DATASIZE)
	print("channel low is ",chn_low)
	print("channel high is ",chn_hgh)
	print("no. of channel is ",channel)
	ychannel, _, _, _, _, _ = read_file(file_in[0], PATH, DATASIZE, chn_low, chn_hgh, channel)
	print("out no. of channel is ",ychannel)
	##########################################################################

	final_res  = np.zeros((nfiles,len(ychannel)))
	final_dat  = np.zeros((nfiles,len(ychannel)))
	final_fit  = np.zeros((nfiles,len(ychannel)))
	final_flg  = np.zeros((nfiles,len(ychannel)))
	final_err  = np.zeros((nfiles,len(ychannel)))

	##########################################################################

	count = 0

	for nfiles_i in range(0, nfiles):

		print ("\nProcessing Night " + str(count+1) + " of " + str(nfiles) + "...")

		##################################################################
		#
		# processing cross (12) data first

		ychannel, ydata, yflag, yres, yfit, yerror = \
				read_file(file_in[nfiles_i], PATH, DATASIZE, chn_low, chn_hgh, channel)
		
		print("file being read is: ",file_in[nfiles_i])
		
		print("yfit is",yfit)
		print("yfit shape is",yfit.shape)
		print("number of non zero values are",np.count_nonzero(yfit))		

		print("yflag is",yflag)
		print("yflag shape is",yflag.shape)
		print("number of non zero values are",np.count_nonzero(yflag))		
		
		
		print("ydata is",ydata)
		print("ydata shape is",ydata.shape)
		print("number of non zero values are",np.count_nonzero(ydata))	
		
		print("yres is",yres)
		print("yres shape is",yres.shape)
		print("number of non zero values are",np.count_nonzero(yres))	
		
		print("ychannel is",ychannel)
		print("ychannel shape is",ychannel.shape)
		print("number of non zero values are",np.count_nonzero(ychannel))
		
		print("yerror is",yerror)
		print("yerror shape is",yerror.shape)
		print("number of non zero values are",np.count_nonzero(yerror))
		
		

		##################################################################
	
		avg_res  = np.zeros(len(ychannel))
		avg_dat  = np.zeros(len(ychannel))
		avg_fit  = np.zeros(len(ychannel))
		avg_err  = np.zeros(len(ychannel))
		flg      = np.ones (len(ychannel)) #By default all channels are good

		##################################################################

		bad_spec = []
		for j in range(0, len(ydata)):
			print("number of non zero value in yflag", np.count_nonzero(yflag[j]))
			if np.sum(yflag[j])==0:
				bad_spec.append(j)
		print("shape of yflag is",np.shape(yflag))
		good_spec = np.delete(range(0, len(ydata)), bad_spec)
		print("good_spec is",good_spec)
		print("good_spec shape is",np.shape(good_spec))
		print("number of non zero values are",np.count_nonzero(good_spec))

		##################################################################

		for i in range(0, len(ychannel)): 

			# With weights:	
			ind = (np.where((yflag[:,i]==1) & (yfit[:,i]!=0.0)))[0]

			if len(ind)==0:   	
				avg_dat[i]  = 0
				avg_res[i]  = 0
				avg_err[i]  = 0
				flg[i]      = 0
				avg_fit[i]  = np.sum(yfit[good_spec, i]) /(len(good_spec))
				continue
								
			avg_res[i]  = np.sum(1.0/yerror[ind, i]**2 * yres[ind, i])/(np.sum(1.0/yerror[ind, i]**2))
			avg_dat[i]  = np.sum(1.0/yerror[ind, i]**2 * ydata[ind, i])/(np.sum(1.0/yerror[ind, i]**2))
			avg_fit[i]  = np.sum(yfit[good_spec, i]) /(len(good_spec))
			avg_err[i]  = 1.0/np.sqrt(np.sum(1.0/yerror[ind, i]**2))

		##################################################################

		err          	       = avg_err
		dat          	       = avg_fit + avg_res
		frq	     	       	   = c_f(ychannel)

		yres         	       = np.zeros(len(ychannel))
		yfit         	       = np.zeros(len(ychannel))
		f_flag         	       = np.zeros(len(ychannel))

		final_res[nfiles_i]    = avg_res	
		final_fit[nfiles_i]    = avg_fit	
		final_dat[nfiles_i]    = dat	
		final_flg[nfiles_i]    = flg
		final_err[nfiles_i]    = err
	
		##################################################################

		count 		       = count + 1

		##################################################################

	##########################################################################

	final_res_c = np.zeros((nfiles,(DATASIZE)), dtype=np.complex64)
	final_dat_c = np.zeros((nfiles,(DATASIZE)), dtype=np.complex64)
	final_fit_c = np.zeros((nfiles,(DATASIZE)), dtype=np.complex64)
	final_flg_c = np.zeros((nfiles,(DATASIZE)), dtype=np.complex64)
	final_err_c = np.zeros((nfiles,(DATASIZE)), dtype=np.complex64)

	##########################################################################

	for i in range(0, len(channel)):
		for j in range(0, len(ychannel)):
			if channel[i]==ychannel[j]: #Assuming same frequency range for all files
				final_res_c[:,i] = final_res[:,j] 
				final_fit_c[:,i] = final_fit[:,j] 	
				final_dat_c[:,i] = final_dat[:,j]    
				final_flg_c[:,i] = final_flg[:,j] 
				final_err_c[:,i] = final_err[:,j] 

    
    
	savetxt('data_obs.csv',final_dat_c.flatten(),delimiter=',')
	k=np.arange(0,8193)
	f=k*500/16384;
	plt.figure()
	plt.xlim([50,100])
	plt.ylim([120,131])
	plt.plot(f,abs(final_dat_c.flatten()),'r-')
	plt.xlabel('Freq_in_MHz')
	plt.ylabel('Amplitude')
	plt.title('Data')
	plt.savefig('Data.png')
	plt.show()
	plt.figure()
	plt.xlim([50,100])
	plt.ylim([120,131])
	plt.plot(f,abs(final_fit_c.flatten()),'r-')
	plt.xlabel('Freq_in_MHz')
	plt.ylabel('Amplitude')
	plt.title('Fit')
	plt.savefig('Fit.png')
	plt.show()
	plt.figure()
	plt.xlim([50,100])
	plt.ylim([0,4])
	plt.plot(f,abs(final_res_c.flatten()),'r-')
	plt.xlabel('Freq_in_MHz')
	plt.ylabel('Amplitude')
	plt.title('Residual')
	plt.savefig('Residual.png')
	plt.show()
	
	print("==========")
	print(np.shape(final_dat_c))
	print(type(final_dat_c))
	print(np.count_nonzero(final_dat_c))
	print("==========")
	##########################################################################

	data_res  = np.ma.array(final_res_c, mask=np.logical_not(final_flg_c), dtype=np.complex64)
	data_dat  = np.ma.array(final_dat_c, mask=np.logical_not(final_flg_c), dtype=np.complex64)
	data_fit  = np.ma.array(final_fit_c, mask=np.logical_not(final_flg_c), dtype=np.complex64)
	data_err  = np.ma.array(final_err_c, mask=np.logical_not(final_flg_c), dtype=np.complex64)
	savetxt('data_residual.csv',final_res_c.flatten(),delimiter=',')
	savetxt('data_fit.csv',final_fit_c.flatten(),delimiter=',')
    

	##########################################################################
		
	uvi       = a.miriad.UV(PATH + file_in[0])
	uv        = a.miriad.UV(PATH + file_out, 'new')
	uv.init_from_uv(uvi)
	del(uvi)

    

	##########################################################################

	for i in range(0,nfiles):
		uvi = a.miriad.UV(PATH + file_in[i])
		for preamble, data in uvi.all():
			uv.copyvr(uvi)
			uv['source'] = uvi['source']
			p_source = uvi['source']
			print ("Write record for source : ",p_source)
			uv['pol']= -1  # Data in RR polarization
			uv.write(preamble, data_dat[i])
			uv['pol']= -2  # Fit in LL polarization 
			uv.write(preamble, data_fit[i])
			uv['pol']= -3  # Residuals in RL polarization
			uv.write(preamble, data_res[i])
			print(ma.getdata(data_dat[i]))
			uv['pol']= -4  # Errors in LR polarization
			uv.write(preamble, data_err[i])
			break
		del(uvi)
    
	del(uv)
	
	##########################################################################
	print("\nTotal time: --- %s seconds ---" % (time.time() - start_time_beg))
	##########################################################################

if __name__ == "__main__":
	main(sys.argv[1:])

