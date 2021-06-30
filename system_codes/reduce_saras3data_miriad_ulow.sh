
#!/bin/bash -f
TIMESTAMP=`date +%Y-%m-%d_%H-%M-%S`
echo $TIMESTAMP

for  srcname in  2020-10-20_092713
do
  s3pprs  "vis = ${srcname}.dat" "out = ${srcname}.pprs" \
 			"tol = 2" "callst = 1" "sitelatitude= 13.992163" "sitelongitude= 74.876027"
  
  s3calib "vis = ${srcname}.pprs" "out = ${srcname}.calib" "phasesw = -1" "tcal = 630"

  s3flag "vis = ${srcname}.calib" \
		   "out = ${srcname}.s3flag" \
		   "chnl = 1312" \
		   "chnh = 2950" \
		   "chnm = 2132" \
		   "polyorder = 10" \
		   "htoll1 = 3" \
		   "htoll2 = 2" \
		   "maxstd = 18" \
 		   "autoh = 5000" \
 		   "autol = 100" \
		   "crossh = 100" \
		   "crossl = -100"
TIMESTAMP=`date +%Y-%m-%d_%H-%M-%S`
echo $TIMESTAMP
 #	python3.7 /home/pi/Desktop/PRATUSH-master/system_codes/s3pyflag_ulow.py -i ${srcname}.s3flag

 # 	python3.7 /home/pi/Desktop/PRATUSH-master/system_codes/s3pymedian_v7.py -i ${srcname}.s3flagpy 
 #	python3.7  /home/pi/Desktop/PRATUSH-master/system_codes/s3pymedian_v8.py -i ${srcname}.s3flagpy -l ${srcname}CH_10Mar_21Mar.txt

 #	python3.7  /home/pi/Desktop/PRATUSH-master/system_codes/s3RFIflag_v3.py -i ${srcname}.s3flagpy.medianpyL -o CH

### Then do a UVCAT to create a single file with extension .RFIflag

	#python3.7  /home/pi/Desktop/PRATUSH-master/system_codes/AvgRecFiles_v7.py -i ${srcname}.s3flagpy.medianpyL.RFIflag

 python3.7 /home/pi/Desktop/PRATUSH-master/system_codes/s3pyflag_ulow.py -i 2020-10-20_092713.s3flag
 python3.7 /home/pi/Desktop/PRATUSH-master/system_codes/s3pymedian_v7.py -i 2020-10-20_092713.s3flagpy -l CH_10Mar_21Mar.txt
 python3.7  /home/pi/Desktop/PRATUSH-master/system_codes/s3RFIflag_v3.py -i 2020-10-20_092713.s3flagpy.medianpyL -o CH
 python3.7  /home/pi/Desktop/PRATUSH-master/system_codes/AvgRecFiles_v7.py -i 2020-10-20_092713.s3flagpy.medianpyL.RFIflag



done

