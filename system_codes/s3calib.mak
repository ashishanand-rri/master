s3calib: s3calib_v3.o 
	 cc -o s3calib s3calib_v3.o \
	 -L/usr/lib \
	 -L/home/pi/Desktop/miriad/linux/lib \
	 -lmir \
	 -lm
s3calib_v3.o : s3calib_v3.c
	  cc -c -DANSI -I/usr/include -I/home/pi/Desktop/PRATUSH-master/system_codes s3calib_v3.c


