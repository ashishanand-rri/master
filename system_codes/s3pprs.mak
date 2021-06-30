s3pprs: s3pprs_v1.o hampel.o nrutil.o cal_lst_v3.o
	 cc -o s3pprs s3pprs_v1.o hampel.o nrutil.o cal_lst_v3.o \
	 -L/usr/lib \
	 -L/home/pi/Desktop/miriad/linux/lib \
	 -lmir \
	 -lm
s3pprs_v1.o : s3pprs_v1.c
	  cc -c -DANSI -I/usr/include -I/home/pi/Desktop/PRATUSH-master/system_codes s3pprs_v1.c
hampel.o : hampel.c
	  cc -c -DANSI -I/usr/include -I/home/pi/Desktop/PRATUSH-master/system_codes  hampel.c
nrutil.o : nrutil.c
		cc -c -DANSI nrutil.c
cal_lst_v3.o : cal_lst_v3.c
		cc -c -DANSI cal_lst_v3.c
