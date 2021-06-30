import numpy as np
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


print("Inside my code")

def main(argv):


	print("inside code of main")
	file_in=[]	
	file_out = ''
	try:
		opts, args = getopt.getopt(argv,"hi:o:","ifile=")
	except getopt.GetoptError:
		print ('avg_output_reader.py -i <inputfile1>')
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print ('avg_output_reader.py -i <inputfile1> ')
			sys.exit()
		elif opt in ("-i", "--ifile"):
			file_in.append(arg)
			for i in range (len(args)):
				file_in.append(args[i])
	file_out   = file_in[0] + 'avg'

	print  ('\nInput file(s) is : ',file_in)
	print  ('Output file is: ',file_out)

	nfiles = len(file_in)
	
	
	DATASIZE = 8193
	channel       = np.asarray(range(1, DATASIZE+1))  # 1, 2, 3....8194

	chn_low       = 1802 # 1638    
	chn_hgh       = 2788 # 3278 

	
	

	PATH          = '/home/pi/Desktop/PRATUSH-master/system_codes/'	

	for nfiles_i in range(0, nfiles):

		ychannel, ydata, yflag, yres, yfit, yerror = read_file(file_in[nfiles_i], PATH, DATASIZE, chn_low, chn_hgh, channel)
				
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
