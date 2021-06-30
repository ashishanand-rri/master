import numpy as np
import aipy as a
import matplotlib.pyplot as plt
import numpy.ma as ma
import numpy.polynomial.legendre as leg
import sys
import os
import time
import itertools
from scipy.optimize import fmin
from scipy import interpolate
import itertools
import time

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

def read_file(file_in, PATH,  DATASIZE, chn_low, chn_hgh, channel):
	
	print("=========================Inside read_file function of read_miriad.py script===========================")
	uvi = a.miriad.UV(PATH + file_in)
	uvi.select('polarization', -3, 0,  include=True)
	# The first read is just to find the number of records to open appropriately sized arrays
	no_records=0
	for preamble, data in uvi.all():
		no_records+=1  # To get the number of records
	uvi.rewind()

#	print "The input file has ",no_records," records "

	# Initialize 2D arrays to store data, fit and flags
	# These are for writing to MIRIAD

	p_sn=[]
	p_source=[]
	p_data = np.zeros((no_records, DATASIZE), dtype=np.complex64)  # real+imag part
	p_flag = np.zeros((no_records, DATASIZE), dtype=np.int64) # integer flags

	i=0
	for preamble, data in uvi.all():
		p_sn.append(i) # Serial number for each record; starts with 0
		p_source.append(uvi['source']) # Source list
		flags=np.logical_not(data.mask) # Miriad flags: 1 is gooddata, Mask: False is gooddata
		
		p_data[i]=data.data #Data array
		p_flag[i]=flags # Miriad type Flags array
		i += 1
	print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ flags are",p_flag)
	print("dimension of flag is",p_flag.shape)
	print("no. of non zeros in p_flag is ",np.count_nonzero(p_flag))
	#To get the starting point right and set store_start to be the starting index
	for i in p_sn:
		if (p_source[i]=='SPCD' and p_source[i+1]=='SPCD' and p_source[i+2]=='SPCD' \
				and p_source[i+3]=='SPCF'):
			store_start=i 
			break

	print ("Start flagging from rec num. ",store_start," from where we have 3 SPCD recs followed by SPCF")

	# Extract the indices for required frequency range chn_low to chn_hgh using first record
	ydata  = p_data[store_start].real # Real part of cross correlation
	yflag  = p_flag[store_start]      # Flag for each record
	yerror = p_data[store_start+3].imag # STD for each record
	_, _, _, _, ilow, ihigh = select_freq_1d(channel, ydata, yerror, yflag, chn_low, chn_hgh) 

	#print "Will use index range: ",ilow,ihigh," to select the channels: ",chn_low,chn_hgh

	lengthfromstart = len(p_data) - store_start
	residuallength = (len(p_data) - store_start)%4
	SPEC_END   = int((lengthfromstart - residuallength)/4)

	ychannel  = channel[ilow:ihigh+1]
	ydata     = np.zeros((SPEC_END, len(ychannel)), dtype=np.float64)
	yres      = np.zeros((SPEC_END, len(ychannel)), dtype=np.float64)
	yfit      = np.zeros((SPEC_END, len(ychannel)), dtype=np.float64)
	yflag     = np.zeros((SPEC_END, len(ychannel)), dtype=np.int32)
	yerror    = np.zeros((SPEC_END, len(ychannel)), dtype=np.float64)

	#######################################################################
	j = 0
	for i in range(store_start, 4*SPEC_END, 4):
		yres[j]	  = p_data[i][ilow:ihigh+1].real
		yflag[j]  = p_flag[i][ilow:ihigh+1]
		j += 1
	del(uvi)
	
	print("yflag is ",yflag)
	print("%%%%%%%%%%%%%%%%%%%%%%%%%%%% no. of non-zero yflag between ilow and ihigh is",np.count_nonzero(yflag))
	print("/n")
	#######################################################################
	uvr = a.miriad.UV(PATH + file_in)
	uvr.select('polarization', -1, 0,  include=True)
	i=0
	for preamble, data in uvr.all():
		p_data[i]=data.data #Data array
		i += 1

	j = 0
	for i in range(store_start, 4*SPEC_END, 4):
		ydata[j]  = p_data[i][ilow:ihigh+1].real
		yerror[j] = p_data[i][ilow:ihigh+1].imag
		j += 1
	
	del(uvr)
	print("ydata inside read_miriad script is ", ydata)
	#######################################################################
	uvs = a.miriad.UV(PATH + file_in)
	uvs.select('polarization', -2, 0,  include=True)

	i=0
	for preamble, data in uvs.all():
		p_data[i]=data.data #Data array
		i += 1

	j = 0
	for i in range(store_start, 4*SPEC_END, 4):
		yfit[j] = p_data[i][ilow:ihigh+1].real
		j += 1

	del(uvs)
	print("========================= Exiting read_file function of read_miriad.py script===========================")
	#######################################################################
	return ychannel, ydata, yflag, yres, yfit, yerror


def read_file11(file_in, PATH,  DATASIZE, chn_low, chn_hgh, channel):

	uvi = a.miriad.UV(PATH + file_in)
	uvi.select('polarization', -3, 0,  include=True)
	# The first read is just to find the number of records to open appropriately sized arrays
	no_records=0
	for preamble, data in uvi.all():
		no_records+=1  # To get the number of records
	uvi.rewind()

#	print "The input file has ",no_records," records "

	# Initialize 2D arrays to store data, fit and flags
	# These are for writing to MIRIAD

	p_sn=[]
	p_source=[]
	p_data = np.zeros((no_records, DATASIZE), dtype=np.complex64)  # real+imag part
	p_flag = np.zeros((no_records, DATASIZE), dtype=np.int64) # integer flags

	i=0
	for preamble, data in uvi.all():
		p_sn.append(i) # Serial number for each record; starts with 0
		p_source.append(uvi['source']) # Source list
		flags=np.logical_not(data.mask) # Miriad flags: 1 is gooddata, Mask: False is gooddata
		p_data[i]=data.data #Data array
		p_flag[i]=flags # Miriad type Flags array
		i += 1

	#To get the starting point right and set store_start to be the starting index
	for i in p_sn:
		if (p_source[i]=='SPCD' and p_source[i+1]=='SPCD' and p_source[i+2]=='SPCD' \
				and p_source[i+3]=='SPCF'):
			store_start=i 
			break

	#print "Start flagging from rec num. ",store_start," from where we have 3 SPCD recs followed by SPCF"

	# Extract the indices for required frequency range chn_low to chn_hgh using first record
	ydata  = p_data[store_start].real # Real part of cross correlation
	yflag  = p_flag[store_start]      # Flag for each record
	yerror = p_data[store_start+3].imag # STD for each record
	_, _, _, _, ilow, ihigh = select_freq_1d(channel, ydata, yerror, yflag, chn_low, chn_hgh) 

	#print "Will use index range: ",ilow,ihigh," to select the channels: ",chn_low,chn_hgh

	lengthfromstart = len(p_data) - store_start
	residuallength = (len(p_data) - store_start)%4
	SPEC_END   = (lengthfromstart - residuallength)/4 

	ychannel  = channel[ilow:ihigh+1]
	ydata     = np.zeros((SPEC_END, len(ychannel)), dtype=np.float64)
	yres      = np.zeros((SPEC_END, len(ychannel)), dtype=np.float64)
	yfit      = np.zeros((SPEC_END, len(ychannel)), dtype=np.float64)
	yflag     = np.zeros((SPEC_END, len(ychannel)), dtype=np.int32)
	yerror    = np.zeros((SPEC_END, len(ychannel)), dtype=np.float64)

	#######################################################################
	j = 0
	for i in range(store_start, 4*SPEC_END, 4):
		yres[j]	  = p_data[i+1][ilow:ihigh+1].real
		yflag[j]  = p_flag[i+1][ilow:ihigh+1]
		j += 1
	del(uvi)
	#######################################################################
	uvr = a.miriad.UV(PATH + file_in)
	uvr.select('polarization', -1, 0,  include=True)
	i=0
	for preamble, data in uvr.all():
		p_data[i]=data.data #Data array
		i += 1

	j = 0
	for i in range(store_start, 4*SPEC_END, 4):
		ydata[j]  = p_data[i+1][ilow:ihigh+1].real
		yerror[j] = p_data[i+3][ilow:ihigh+1].imag
		j += 1

	del(uvr)
	#######################################################################
	uvs = a.miriad.UV(PATH + file_in)
	uvs.select('polarization', -2, 0,  include=True)

	i=0
	for preamble, data in uvs.all():
		p_data[i]=data.data #Data array
		i += 1

	j = 0
	for i in range(store_start, 4*SPEC_END, 4):
		yfit[j] = p_data[i+1][ilow:ihigh+1].real
		j += 1

	del(uvs)
	#######################################################################
	return ychannel, ydata, yflag, yres, yfit, yerror

def read_file22(file_in, PATH,  DATASIZE, chn_low, chn_hgh, channel):

	uvi = a.miriad.UV(PATH + file_in)
	uvi.select('polarization', -3, 0,  include=True)
	# The first read is just to find the number of records to open appropriately sized arrays
	no_records=0
	for preamble, data in uvi.all():
		no_records+=1  # To get the number of records
	uvi.rewind()

#	print "The input file has ",no_records," records "

	# Initialize 2D arrays to store data, fit and flags
	# These are for writing to MIRIAD

	p_sn=[]
	p_source=[]
	p_data = np.zeros((no_records, DATASIZE), dtype=np.complex64)  # real+imag part
	p_flag = np.zeros((no_records, DATASIZE), dtype=np.int64) # integer flags

	i=0
	for preamble, data in uvi.all():
		p_sn.append(i) # Serial number for each record; starts with 0
		p_source.append(uvi['source']) # Source list
		flags=np.logical_not(data.mask) # Miriad flags: 1 is gooddata, Mask: False is gooddata
		p_data[i]=data.data #Data array
		p_flag[i]=flags # Miriad type Flags array
		i += 1

	#To get the starting point right and set store_start to be the starting index
	for i in p_sn:
		if (p_source[i]=='SPCD' and p_source[i+1]=='SPCD' and p_source[i+2]=='SPCD' \
				and p_source[i+3]=='SPCF'):
			store_start=i 
			break

	print ("Start flagging from rec num. ",store_start," from where we have 3 SPCD recs followed by SPCF")

	# Extract the indices for required frequency range chn_low to chn_hgh using first record
	ydata  = p_data[store_start].real # Real part of cross correlation
	yflag  = p_flag[store_start]      # Flag for each record
	yerror = p_data[store_start+3].imag # STD for each record
	_, _, _, _, ilow, ihigh = select_freq_1d(channel, ydata, yerror, yflag, chn_low, chn_hgh) 

	#print "Will use index range: ",ilow,ihigh," to select the channels: ",chn_low,chn_hgh

	lengthfromstart = len(p_data) - store_start
	residuallength = (len(p_data) - store_start)%4
	SPEC_END   = (lengthfromstart - residuallength)/4 

	ychannel  = channel[ilow:ihigh+1]
	ydata     = np.zeros((SPEC_END, len(ychannel)), dtype=np.float64)
	yres      = np.zeros((SPEC_END, len(ychannel)), dtype=np.float64)
	yfit      = np.zeros((SPEC_END, len(ychannel)), dtype=np.float64)
	yflag     = np.zeros((SPEC_END, len(ychannel)), dtype=np.int32)
	yerror    = np.zeros((SPEC_END, len(ychannel)), dtype=np.float64)

	#######################################################################
	j = 0
	for i in range(store_start, 4*SPEC_END, 4):
		yres[j]	  = p_data[i+2][ilow:ihigh+1].real
		yflag[j]  = p_flag[i+2][ilow:ihigh+1]
		j += 1
	del(uvi)
	#######################################################################
	uvr = a.miriad.UV(PATH + file_in)
	uvr.select('polarization', -1, 0,  include=True)
	i=0
	for preamble, data in uvr.all():
		p_data[i]=data.data #Data array
		i += 1

	j = 0
	for i in range(store_start, 4*SPEC_END, 4):
		ydata[j]  = p_data[i+2][ilow:ihigh+1].real
		yerror[j] = p_data[i+3][ilow:ihigh+1].imag
		j += 1

	del(uvr)
	#######################################################################
	uvs = a.miriad.UV(PATH + file_in)
	uvs.select('polarization', -2, 0,  include=True)

	i=0
	for preamble, data in uvs.all():
		p_data[i]=data.data #Data array
		i += 1

	j = 0
	for i in range(store_start, 4*SPEC_END, 4):
		yfit[j] = p_data[i+2][ilow:ihigh+1].real
		j += 1

	del(uvs)
	#######################################################################
	return ychannel, ydata, yflag, yres, yfit, yerror
