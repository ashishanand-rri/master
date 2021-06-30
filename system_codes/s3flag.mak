s3flag: s3flag.o nrutil.o amoeba.o amotry.o select_nr.o
	 cc -o s3flag s3flag.o nrutil.o amoeba.o amotry.o select_nr.o \
	 -L/usr/lib \
	 -L/home/pi/Desktop/miriad/linux/lib \
	 -lmir \
	 -lm
s3flag.o : s3flag_v3.c
	cc -c -o s3flag.o -DANSI -I/usr/include -I/home/pi/Desktop/PRATUSH-master/system_codes s3flag_v3.c
nrutil.o : nrutil.c
	cc -c -DANSI nrutil.c
amotry.o : amotry.c
	cc -c -DANSI -I/home/pi/Desktop/PRATUSH-master/system_codes  amotry.c
amoeba.o : amoeba.c
	cc -c -DANSI -I/home/pi/Desktop/PRATUSH-master/system_codes amoeba.c
select_nr.o : select_nr.c
	cc -c -DANSI -I/home/pi/Desktop/PRATUSH-master/system_codes select_nr.c
