/*= S3CALIB - Process the saras3 data file after averaging within scans */
/*& Ravi */
/*: Utility */
/*+
 S3CALIB* reads the miriad uv dataset containing cal and obs scans.
 Does a simple calibration of the cross spectra 
 assuming no interference.  Writes out a calibrated file.

 v2 is for 6-state optical front end with codes cycling through
 OBS00 OBS11 CAL00 CAL01 CAL10 CAL11

 v3 calibrates the difference between (OBS00-OBS11) and (CAL00-CAL01)
 to remove additives from the reference port.  A double difference!
*/
/*@ vis */
/*         Name of the uvdataset; no defaults */
/*@ out */
/*         Name of the calibrated uvdataset; no defaults */
/*@ tcal */
/*         Tcal in K; default 446.0 */
/*@ phasesw */
/*         Switch phase on CALA0, give (+1/-1); default -1 */
/*-- */
/*
c***************************************************************
c
c History
c
c Ravi 14Aug17 original version from s2calib_v7.c
c
c***************************************************************
*/

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <maxdimc.h>
#include <miriad.h>

#define PI 3.14159265
#define DATASIZE 8193
#define INT_SPEC 67.108864  // ms of integration time per record of raw dat file


float timedifference_msec(struct timeval t0,struct timeval t1)
{
	return(t1.tv_sec-t0.tv_sec)*1000.0f + (t1.tv_usec-t0.tv_usec)/1000.0f;
}

int main(argc,argv)
int argc;
char *argv[];
{
	
	
struct timeval t_zero;
struct timeval t_one;
float elapsed;
gettimeofday(&t_zero,0);	
	
const char version[30]="s3calib: version 3 06mar18";
char visfile[50],outfile[50],text[50];
char status[30]="old";
int i,j,k,nn,ind;
int tno,tno_in;
int ndata,nread;
int iwin,datasize;
float a_12,b_12,c_12,d_12,e_12,f_12,ai_12,bi_12,ci_12;
float r_12,ri_12;
float di_12,ei_12,fi_12,g_12,gi_12,h_12,hi_12,cc_12,cci_12;
float a_11,b_11,c_11,d_11,e_11,f_11,g_11,gi_11,h_11,hi_11,cc_11,cci_11;
float a_22,b_22,c_22,d_22,e_22,f_22,g_22,gi_22,h_22,hi_22,cc_22,cci_22;
float aa_12,aai_12,bb_12,bbi_12,gg_12,ggi_12,hh_12,hhi_12;
float aa_11,aai_11,bb_11,bbi_11,gg_11,ggi_11,hh_11,hhi_11;
float aa_22,aai_22,bb_22,bbi_22,gg_22,ggi_22,hh_22,hhi_22;
float bb_mag_12, aa_mag_12, spect_array_mag, t1, t2, t3;
double del_ref00,del_ref11,del_cal00,del_cal11;
double del_obsMref, del_calMref, del_ref, mag_obsMref, mag_calMref;


unsigned long ncal00,ncal01,nobs0,ncal10,ncal11,nobs1,nspect,ispect, nflag;
unsigned long nlst,ilst;
int phasesw;

float tcal;
float cal01_11_array[2*DATASIZE],obs1_11_array[2*DATASIZE];
float cal00_11_array[2*DATASIZE],obs0_11_array[2*DATASIZE];
float cal10_11_array[2*DATASIZE],cal11_11_array[2*DATASIZE];
float cal01_12_array[2*DATASIZE],obs1_12_array[2*DATASIZE];
float cal00_12_array[2*DATASIZE],obs0_12_array[2*DATASIZE];
float cal10_12_array[2*DATASIZE],cal11_12_array[2*DATASIZE];
float cal01_22_array[2*DATASIZE],obs1_22_array[2*DATASIZE];
float cal00_22_array[2*DATASIZE],obs0_22_array[2*DATASIZE];
float cal10_22_array[2*DATASIZE],cal11_22_array[2*DATASIZE];
    
float spect_array0[2*DATASIZE],spect_array1[2*DATASIZE],spect_array_cal[2*DATASIZE];
float spect_array_12[2*DATASIZE];
float spect_array_11[2*DATASIZE];
float spect_array_22[2*DATASIZE];
int flags_array[6][DATASIZE];
double flags_store[6][2*DATASIZE];
double lst_store[6];

// Miriad-related variables
int flags[DATASIZE];
float flags_final[2*DATASIZE];
float data[2*DATASIZE];
double preamble_11,preamble_12,preamble_22,preamble[4];
char source_name[30];
char velocity_type[9] = {'V','E','L','O','-','O','B','S','\0'};
long int ispec;
int fname_length;
int cnumber;
int var_ivalue[1];
static int npoints;
int antenna1,antenna2,antenna3;
int access_result;
float baseline_value[1];
double p1,p2,p3;
double ant_az[3],ant_el[3];
double freq_channel1[1];
double freq_inc[1];
double time_var[1];
double coord_var[2];
double site_latitude[1],site_longitude[1],longitude_radians;
float var_veldop[1];
float var_vsource[1];
double var_restfreq[1];
double sra[1],sdec[1];
double lo1[1],lo2[1],freq[1],freqif[1];
int mount[1];
float evector[1];
int nnpol,npol[1];
float jyperk[1],inttime[1],epoch[1];
double antpos[9];
float tpower[1];
float wfreq[1],wwidth[1];
int time_info[6][16];
int time_info1[6];
long int jjdd,yyyy,mmmm,dddd;
double julian_date;
double lst;
float old_time,old_day,old_year;
float current_time,current_day,current_year;

char pre_source[30];
int warn_flag;
int ant1,ant2,ant_code;

double time_obs00, time_obs11;
double del_obs00, del_obs11, del_cala0, del_cala1, n1a,n1b, n2a, n2b,d1,d2, n2, del_cal, del_obs, del_cal_obs;
double del_obs_sq[DATASIZE], del_cal_sq[DATASIZE];
float P_cal_step;
double P_obs_dif_sq[DATASIZE], del_p_calib_sq[DATASIZE];


/* key calls to get user inputs */
keyini_c(argc,argv);
keya_c("vis",visfile,"");
keya_c("out",outfile,"");
keyr_c("tcal",&tcal,446.00);
keyi_c("phasesw",&phasesw,-1);
keyfin_c();
	
puts(version);

if(phasesw != 1 && phasesw != -1) 
	{
	printf("error in phasesw input\n"); 
	return(0);
	}

uvopen_c(&tno_in,visfile,status);
uvopen_c(&tno,outfile,"new");
hisopen_c(tno,"write");

datasize=DATASIZE;

/* setup the output vis file */
wrhda_c(tno,"obstype","crosscorrelation");
uvputvra_c(tno,"source","zenith sky");
uvputvra_c(tno,"operator","eor");
uvputvra_c(tno,"version",version);
sra[0]=0.0;
sdec[0]=0.0;
uvputvrd_c(tno,"ra",sra,1);
uvputvrd_c(tno,"obsra",sra,1);
uvputvrd_c(tno,"dec",sdec,1);
uvputvrd_c(tno,"obsdec",sdec,1);
lo1[0]=0.250;
lo2[0]=0.0;
freq[0]=0.125;
freqif[0]=0.0;
uvputvrd_c(tno,"lo1",lo1,1);
uvputvrd_c(tno,"lo2",lo2,1);
uvputvrd_c(tno,"freq",freq,1);
uvputvrd_c(tno,"freqif",freqif,1);
mount[0]=0;
evector[0]=0.0;
uvputvri_c(tno,"mount", mount, 1);
uvputvrr_c(tno,"evector",evector,1);
uvputvrr_c(tno,"chi",evector,1);
uvputvra_c(tno,"telescop","eor");
jyperk[0]=1.0;
uvputvrr_c(tno,"jyperk",jyperk,1);
inttime[0]=1.0;
uvputvrr_c(tno,"inttime",inttime,1);
epoch[0]=2000.0;
uvputvrr_c(tno,"epoch",epoch,1);
nnpol=1;
wrhdi_c(tno,"npol",nnpol);
antpos[0]=0.0;
antpos[1]=1.0;
antpos[2]=2.0;
antpos[3]=0.0;
antpos[4]=0.0;
antpos[5]=0.0;
antpos[6]=0.0;
antpos[7]=0.0;
antpos[8]=0.0;
uvputvrd_c(tno,"antpos",antpos,9);
uvputvrd_c(tno,"time",time_var,1);
uvputvrd_c(tno,"ut",time_var,1);
nn=0;
p1=0.0;
p2=0.0;
p3=0.0;
uvset_c(tno,"corr","r",nn,p1,p2,p3);  
/* write floats and not scaled integers */
var_ivalue[0]=3;
uvputvri_c(tno,"nants",var_ivalue,1);
var_ivalue[0]=8193;
uvputvri_c(tno,"nchan",var_ivalue,1);
var_ivalue[0]=1;
uvputvri_c(tno,"npol",var_ivalue,1);
var_ivalue[0]=1;
uvputvri_c(tno,"nspect",var_ivalue,1);
var_ivalue[0]=-1;
uvputvri_c(tno,"pol",var_ivalue,1);
var_ivalue[0]=8193;
uvputvri_c(tno,"nschan",var_ivalue,1);
var_ivalue[0]=1;
uvputvri_c(tno,"ischan",var_ivalue,1);
var_ivalue[0]=1;
uvputvri_c(tno,"ntpower",var_ivalue,1);
var_ivalue[0]=0;
uvputvri_c(tno,"nwide",var_ivalue,1);
tpower[0]=1.0;
uvputvrr_c(tno,"tpower",tpower,1);
ant_az[0]=0.0;
ant_az[1]=0.0;
ant_az[2]=0.0;
ant_el[0]=90.0;
ant_el[1]=90.0;
ant_el[2]=90.0;
uvputvrd_c(tno,"antaz",ant_az,3);
uvputvrd_c(tno,"antel",ant_el,3);
var_veldop[0]=0.0;
uvputvrr_c(tno,"veldop",var_veldop,1);
var_vsource[0]=0.0;
uvputvrr_c(tno,"vsource",var_vsource,1);
var_restfreq[0]=0.0;
uvputvrd_c(tno,"restfreq",var_restfreq,1);
freq_channel1[0]=0.0; /* GHz */
uvputvrd_c(tno,"sfreq",freq_channel1,1);
freq_inc[0]= 3.0517578125e-05; /* GHz */
uvputvrd_c(tno,"sdf",freq_inc,1);
site_latitude[0] = (double)(13.01333 * PI/180.0);
uvputvrd_c(tno,"latitud",site_latitude,1);
site_longitude[0] = (double)(77.580833 * PI/180.0);
uvputvrd_c(tno,"longitu",site_longitude,1);
antenna1=1;
antenna2=2;
baseline_value[0]=(float)(256*antenna1+antenna2);
uvputvrr_c(tno,"baseline",baseline_value,1);
coord_var[0]=0.0;
coord_var[1]=0.0;
uvputvrd_c(tno,"coord",coord_var,2);
uvputvr_c(tno,1,"veltype",velocity_type,8);	

ndata=DATASIZE;
nn=30;
ncal00=0;
ncal01=0;
nobs0=0;
ncal10=0;
ncal11=0;
nobs1=0;
nspect=0;
nlst=0;
nflag=0;

for(;;)  
    {								// outer loop reading records of file
	uvread_c(tno_in,preamble,data,flags,ndata,&nread);

	if(nread == 0) break;

	uvgetvra_c(tno_in,"source",source_name,nn);
	uvgetvrd_c(tno_in,"lst",time_var,1);

	/* preamble[3] is (double)(256*ant1+ant2) */
	ant_code=(int)preamble[3];
	ant2=ant_code%256;
	ant1=(ant_code-ant2)/256;

	if (ant1==1 && ant2==1)
		{
		preamble_11 = preamble[3];
		if (source_name[0] == 'C' && source_name[3] == '0' && source_name[4] == '0')
			{
			for(i=0;i<2*DATASIZE;++i)
			cal00_11_array[i]=data[i];
			for(j=0; j<30; ++j) 
			pre_source[j] = source_name[j];
			++ncal00;
			++nspect;
			}
		if (source_name[0] == 'C' && source_name[3] == '0' && source_name[4] == '1')
			{			
		    for(i=0;i<2*DATASIZE;++i)
			cal01_11_array[i]=data[i];
		    for(j=0; j<30; ++j) 
			pre_source[j] = source_name[j];
		    ++ncal01;
		    ++nspect;
			}
		if (source_name[0] == 'C' && source_name[3] == '1' && source_name[4] == '0')
			{
			for(i=0;i<2*DATASIZE;++i)
			cal10_11_array[i]=data[i];
			for(j=0; j<30; ++j) 
			pre_source[j] = source_name[j];
			++ncal10;
			++nspect;
			}
		if (source_name[0] == 'C' && source_name[3] == '1' && source_name[4] == '1')
			{			
		    for(i=0;i<2*DATASIZE;++i)
			cal11_11_array[i]=data[i];
		    for(j=0; j<30; ++j) 
			pre_source[j] = source_name[j];
		    ++ncal11;
		    ++nspect;
			}
		if (source_name[0] == 'O' && source_name[4] == '0')
			{	
		    for(i=0;i<2*DATASIZE;++i)
			obs0_11_array[i]=data[i];
		    for(j=0; j<30; ++j) 
			pre_source[j] = source_name[j];
		    ++nobs0;
		    ++nspect;
			}
		if (source_name[0] == 'O' && source_name[4] == '1')
			{	
		    for(i=0;i<2*DATASIZE;++i)
			obs1_11_array[i]=data[i];
		    for(j=0; j<30; ++j) 
			pre_source[j] = source_name[j];
		    ++nobs1;
		    ++nspect;
			}
		}  // end of processing of record with 11

	if (ant1==2 && ant2==2)
		{
		preamble_22 = preamble[3];
		if (source_name[0] == 'C' && source_name[3] == '0' && source_name[4] == '0')
			{	
		    for(i=0;i<2*DATASIZE;++i)
			cal00_22_array[i]=data[i];
		    for(j=0; j<30; ++j) 
			pre_source[j] = source_name[j];
		    ++ncal00;
		    ++nspect;
			}
		if (source_name[0] == 'C' && source_name[3] == '0' && source_name[4] == '1')
			{	
		    for(i=0;i<2*DATASIZE;++i)
			cal01_22_array[i]=data[i];
		    for(j=0; j<30; ++j) 
			pre_source[j] = source_name[j];
		    ++ncal01;
		    ++nspect;
			}
		if (source_name[0] == 'C' && source_name[3] == '1' && source_name[4] == '0')
			{	
		    for(i=0;i<2*DATASIZE;++i)
			cal10_22_array[i]=data[i];
		    for(j=0; j<30; ++j) 
			pre_source[j] = source_name[j];
		    ++ncal10;
		    ++nspect;
			}
		if (source_name[0] == 'C' && source_name[3] == '1' && source_name[4] == '1')
			{	
		    for(i=0;i<2*DATASIZE;++i)
			cal11_22_array[i]=data[i];
		    for(j=0; j<30; ++j) 
			pre_source[j] = source_name[j];
		    ++ncal11;
		    ++nspect;
			}
		if (source_name[0] == 'O' && source_name[4] == '0')
			{
		    for(i=0;i<2*DATASIZE;++i)
			obs0_22_array[i]=data[i];
		    for(j=0; j<30; ++j) 
			pre_source[j] = source_name[j];
		    ++nobs0;
		    ++nspect;
			}
		if (source_name[0] == 'O' && source_name[4] == '1')
			{		
		    for(i=0;i<2*DATASIZE;++i)
			obs1_22_array[i]=data[i];
		    for(j=0; j<30; ++j) 
			pre_source[j] = source_name[j];
		    ++nobs1;
		    ++nspect;
			}
		}  // end of processing or record with 22

	/* Apart from cross correlation, FLAG source also has the same preamble as 12
	 Further, FLAG source will always be AFTER one of the standard four sources
	 Hence pre_source[] is never NULL
	 FLAG source is same for all three antenna combinations (1-1, 1-2, 2-2)
	 Mapping of flags_array and flag_store first index :
	 OBS00 - 0
	 OBS11 - 1
	 CAL00 - 2
	 CAL01 - 3
	 CAL10 - 4
	 CAL11 - 5	
	*/

	if (ant1==1 && ant2==2) 
		{
		preamble_12 = preamble[3];
		if(source_name[0] == 'C' && source_name[3] == '0' && source_name[4] == '0')
			{
			for(i=0;i<2*DATASIZE;++i)
				cal00_12_array[i]=data[i];
			for(j=0; j<30; ++j) 
				pre_source[j] = source_name[j];
			lst_store[nlst] = time_var[0];               // get LST only from cross 
			for (k=0;k<ndata;++k) flags_array[2][k]=flags[k];
			++nlst;
			++ncal00;
			++nspect;
			}
		if(source_name[0] == 'C' && source_name[3] == '0' && source_name[4] == '1')
			{
			for(i=0;i<2*DATASIZE;++i)
				cal01_12_array[i]=data[i];
			for(j=0; j<30; ++j) 
				pre_source[j] = source_name[j];
			lst_store[nlst] = time_var[0]; 
			for (k=0;k<ndata;++k) flags_array[3][k]=flags[k];
			++nlst; 
			++ncal01;
			++nspect;
			}
		if(source_name[0] == 'C' && source_name[3] == '1' && source_name[4] == '0')
			{
			for(i=0;i<2*DATASIZE;++i)
				cal10_12_array[i]=data[i];
			for(j=0; j<30; ++j) 
				pre_source[j] = source_name[j];
			lst_store[nlst] = time_var[0];               // get LST only from cross 
			for (k=0;k<ndata;++k) flags_array[4][k]=flags[k];
			++nlst;
			++ncal10;
			++nspect;
			}
		if(source_name[0] == 'C' && source_name[3] == '1' && source_name[4] == '1')
			{
			for(i=0;i<2*DATASIZE;++i)
				cal11_12_array[i]=data[i];
			for(j=0; j<30; ++j) 
				pre_source[j] = source_name[j];
			lst_store[nlst] = time_var[0]; 
			for (k=0;k<ndata;++k) flags_array[5][k]=flags[k];
			++nlst; 
			++ncal11;
			++nspect;
			}
		if(source_name[0] == 'O' && source_name[4] == '0')
			{
			for(i=0;i<2*DATASIZE;++i)
				obs0_12_array[i]=data[i];
			for(j=0; j<30; ++j) 
				pre_source[j] = source_name[j];
			lst_store[nlst] = time_var[0];
			for (k=0;k<ndata;++k) flags_array[0][k]=flags[k];
			++nlst;
			++nobs0;
			++nspect;
			}
		if(source_name[0] == 'O' && source_name[4] == '1')
			{
			for(i=0;i<2*DATASIZE;++i)
				obs1_12_array[i]=data[i];
			for(j=0; j<30; ++j) 
				pre_source[j] = source_name[j];
			lst_store[nlst] = time_var[0];
			for (k=0;k<ndata;++k) flags_array[1][k]=flags[k];
			++nlst;
			++nobs1;
			++nspect;
			}

		if(source_name[2] == 'A' && source_name[3] == 'G')
			{
			if(pre_source[0] == 'O' && pre_source[4] == '0')
				{
				for(i=0;i<2*DATASIZE;++i) 
					flags_store[0][i]=data[i];
				++nflag;
				++nspect;
				}
			if(pre_source[0] == 'O' && pre_source[4] == '1')
				{
				for(i=0;i<2*DATASIZE;++i) 
					flags_store[1][i]=data[i];
				++nflag;
				++nspect;
				}
			if(pre_source[0] == 'C' && pre_source[3] == '0' && pre_source[4] == '0')
				{
				for(i=0;i<2*DATASIZE;++i) 
					flags_store[2][i]=data[i];
				++nflag;
				++nspect;
				}
			if(pre_source[0] == 'C' && pre_source[3] == '0' && pre_source[4] == '1')
				{
				for(i=0;i<2*DATASIZE;++i) 
					flags_store[3][i]=data[i];
				++nflag;
				++nspect;
				}
			if(pre_source[0] == 'C' && pre_source[3] == '1' && pre_source[4] == '0')
				{
				for(i=0;i<2*DATASIZE;++i) 
					flags_store[4][i]=data[i];
				++nflag;
				++nspect;
				}
			if(pre_source[0] == 'C' && pre_source[3] == '1' && pre_source[4] == '1')
				{
				for(i=0;i<2*DATASIZE;++i) 
					flags_store[5][i]=data[i];
				++nflag;
				++nspect;
				}
			}
		}  // end of processing of record with 12

	if(nspect==24)  // come here if 1*1, 2*2 1*2 "baselines" are read in 
					// with obs00,obs11,cal00,cal01,cal10,cal11 and total 6 flag records 
	{
	if(ncal00!=3 || ncal01!=3 || ncal10!=3 || ncal11!=3 ||nobs0!=3 || nobs1!=3 || nflag!=6 || nlst!=6 ) 
		{ printf("Inconsistent data....exiting \n"); break;}

	/* compute the mean LST */
	time_var[0]=0.0;
	for(ilst=0;ilst<6;++ilst)
	{
	time_var[0] += lst_store[ilst];
	}
	time_var[0] /= (double)6.0;


	/* set the output flag to 1 if all six obs00,obs11,cal00,cal01,cal10,cal11 are good */
	for (k=0;k<ndata;++k)
		{
		if(flags_array[0][k]==0 || flags_array[1][k]==0 || 
				flags_array[2][k]==0 || flags_array[3][k]==0 ||
				flags_array[4][k]==0 || flags_array[5][k]==0) flags[k]=0;
		else flags[k]=1;
		}

	for(i=0;i<DATASIZE;++i) // For cross 12
	    {
		a_12=cal00_12_array[2*i];			// a, ai = CAL00
		b_12=cal01_12_array[2*i];			// b, bi = CAL01
		ai_12=cal00_12_array[2*i+1];		
		bi_12=cal01_12_array[2*i+1];

		aa_12=cal10_12_array[2*i];			// aa, aai = CAL10
		bb_12=cal11_12_array[2*i];			// bb, bbi = CAL11
		aai_12=cal10_12_array[2*i+1];		
		bbi_12=cal11_12_array[2*i+1];

		c_12=obs0_12_array[2*i];			// c, ci = OBS0
		d_12=obs1_12_array[2*i];			// d, di = OBS1
		ci_12=obs0_12_array[2*i+1];
		di_12=obs1_12_array[2*i+1];

		// subtract cold load from hot load to get the CAL STEP

		e_12=aa_12-a_12; 				// aa, aai = CAL-REF
		ei_12=aai_12-ai_12;
		f_12=bb_12-b_12;
		fi_12=bbi_12-bi_12;

		// difference the values in successive switch positions

		c_12 = c_12-d_12;
		ci_12 = ci_12-di_12;  // OBS0-OBS1

		r_12 = a_12-b_12;
		ri_12 = ai_12-bi_12;  // REF0-REF1

		c_12 -= r_12;
		ci_12 -= ri_12;       // Subtract REF from OBS

		e_12 = e_12-f_12;
		ei_12 = ei_12-fi_12;  // CALON-CALOFF

		// divide OBS/CAL and multiply by tcal

		spect_array_12[2*i]=(tcal)*(c_12*e_12+ci_12*ei_12)/(e_12*e_12+ei_12*ei_12);
		spect_array_12[2*i+1]=(tcal)*(ci_12*e_12-c_12*ei_12)/(e_12*e_12+ei_12*ei_12);
	    }	
	source_name[0]='S';
	source_name[1]='P';
	source_name[2]='C';
	source_name[3]='D';
	source_name[4]='\0';
	    
    preamble[3]=preamble_12;
	uvputvr_c(tno,1,"source",source_name,4);
    uvputvrd_c(tno,"lst",time_var,1);
	uvwrite_c(tno,preamble,spect_array_12,flags,nread);
	    
	for(i=0;i<DATASIZE;++i)	// for auto 11
	    {
		a_11=cal00_11_array[2*i];
		b_11=cal01_11_array[2*i];  // CALOFF

		aa_11=cal10_11_array[2*i];
		bb_11=cal11_11_array[2*i];  // CALON

		c_11=obs0_11_array[2*i];
		d_11=obs1_11_array[2*i];    // OBS

		e_11 = aa_11 - a_11;
		f_11 = bb_11 - b_11;        // CALON-CALOFF 

		spect_array_11[2*i]=(tcal)*(c_11+d_11)/(e_11+f_11);
		spect_array_11[2*i+1]=0.0; 
	    }	  
	source_name[0]='S';
	source_name[1]='P';
	source_name[2]='C';
	source_name[3]='D';
	source_name[4]='\0';
	preamble[3]=preamble_11;
	uvputvr_c(tno,1,"source",source_name,4);
	uvputvrd_c(tno,"lst",time_var,1);
	uvwrite_c(tno,preamble,spect_array_11,flags,nread);
  
	for(i=0;i<DATASIZE;++i)	// for auto 22
	    {
		a_22=cal00_22_array[2*i];
		b_22=cal01_22_array[2*i];  // CALOFF

		aa_22=cal10_22_array[2*i];
		bb_22=cal11_22_array[2*i];  // CALON

		c_22=obs0_22_array[2*i];
		d_22=obs1_22_array[2*i];    // OBS

		e_22 = aa_22 - a_22;
		f_22 = bb_22 - b_22;        // CALON-CALOFF 

		spect_array_22[2*i]=(tcal)*(c_22+d_22)/(e_22+f_22);
		spect_array_22[2*i+1]=0.0; 
	    }	  
	source_name[0]='S';
	source_name[1]='P';
	source_name[2]='C';
	source_name[3]='D';
	source_name[4]='\0';
	preamble[3]=preamble_22;
	uvputvr_c(tno,1,"source",source_name,4);
	uvputvrd_c(tno,"lst",time_var,1);
	uvwrite_c(tno,preamble,spect_array_22,flags,nread);
	
	for(i=0;i<DATASIZE;++i)  /* compute the rms in each channel of the calibrated spectra */
		{
		if(flags[i]==0) // set "data" to zero and skip computing rms
			{
			flags_final[2*i]=0;
			flags_final[2*i+1]=0;
			}
		else 			// compute rms for this good channel
			{

	// Continue for a channel only when ALL the 6 values coming from 6 sources are unflagged
	// flags_store can have value from 1-24 below (no. of unflagged channels out of 24)

	// Calculation of equivalent integration time (ms), write as real part of SPCF
	// Equivalent integration time for a source = (flags_store * 67.108864) ms

			time_obs00   =(double) flags_store[0][2*i];
			time_obs11   =(double) flags_store[1][2*i];
			flags_final[2*i]   = (float)(time_obs00 + time_obs11)*INT_SPEC;	  // ms units
						

	// Calculation of variance, imaginary part of SPCF
	// flags_store contains standard deviation for imaginary of cross
	// However, real will also have the same value	

			a_12=cal00_12_array[2*i];			// a, ai = CAL00
			b_12=cal01_12_array[2*i];			// b, bi = CAL01
			aa_12=cal10_12_array[2*i];			// aa, aai = CAL10
			bb_12=cal11_12_array[2*i];			// bb, bbi = CAL11
			ai_12=cal00_12_array[2*i+1];		
			bi_12=cal01_12_array[2*i+1];
			aai_12=cal10_12_array[2*i+1];		
			bbi_12=cal11_12_array[2*i+1];

			c_12=obs0_12_array[2*i];			// c, ci = OBS00
			d_12=obs1_12_array[2*i];			// d, di = OBS11
			ci_12=obs0_12_array[2*i+1];
			di_12=obs1_12_array[2*i+1];

			// subtract cold load from hot load to get the CAL STEP

			e_12=aa_12-a_12; 				// aa, aai = CAL-REF
			ei_12=aai_12-ai_12;
			f_12=bb_12-b_12;
			fi_12=bbi_12-bi_12;

			// difference the values in successive switch positions

			c_12 = c_12-d_12;
			ci_12 = ci_12-di_12;  // OBS0-OBS1

			r_12 = a_12-b_12;
			ri_12 = ai_12-bi_12;  // REF0-REF1

			c_12 -= r_12;
			ci_12 -= ri_12;       // Subtract REF from OBS

			e_12 = e_12-f_12;
			ei_12 = ei_12-fi_12;  // CALON-CALOFF

			del_obs00        = flags_store[0][2*i+1]/sqrt(flags_store[0][2*i]);
			del_obs11        = flags_store[1][2*i+1]/sqrt(flags_store[1][2*i]);
			del_ref00        = flags_store[2][2*i+1]/sqrt(flags_store[2][2*i]);
			del_ref11        = flags_store[3][2*i+1]/sqrt(flags_store[3][2*i]);
			del_cal00        = flags_store[4][2*i+1]/sqrt(flags_store[4][2*i]);
			del_cal11        = flags_store[5][2*i+1]/sqrt(flags_store[5][2*i]);

			del_obsMref = sqrt( pow(del_obs00,2)+pow(del_obs11,2)+pow(del_ref00,2)+pow(del_ref11,2) );
			del_calMref = sqrt( pow(del_cal00,2)+pow(del_cal11,2)+pow(del_ref00,2)+pow(del_ref11,2) );

			del_ref = sqrt( pow(del_ref00,2)+pow(del_ref11,2) );
			mag_obsMref = sqrt( c_12*c_12 + ci_12*ci_12 );
			mag_calMref = sqrt( e_12*e_12 + ei_12*ei_12 );
			spect_array_mag  = sqrt(powf(spect_array_12[2*i],2) + powf(spect_array_12[2*i+1],2));

			t1 = powf(del_obsMref/mag_obsMref, 2);
			t2 = powf(del_calMref/mag_calMref, 2);
			t3 = (del_ref*del_ref/(mag_obsMref*mag_calMref));
			flags_final[2*i+1]=sqrt(t1 + t2 - 2*t3) * spect_array_mag;

			} // done with computing rms for this channel

		} // end of loop through all channels

	    source_name[0]='S';
	    source_name[1]='P';
	    source_name[2]='C';
	    source_name[3]='F';
	    source_name[4]='\0';
	   
        preamble[3]=preamble_12; 	// preamble details, apart from source name, 
									// are same for FLAG and Cross correlation spectra
	    uvputvr_c(tno,1,"source",source_name,4);
        uvputvrd_c(tno,"lst",time_var,1);
	    uvwrite_c(tno,preamble,flags_final,flags,nread);
	   
	    nspect=0;
		nlst=0;
	    ncal00=0;
        ncal01=0;
	    ncal10=0;
        ncal11=0;
	    nobs0=0;
        nobs1=0;
	    nflag=0;

	}  	// written the 11, 22 and 12 SPCD records and an SPCF record 
		// at the end of processing 16 record block


    }  // end of outer loop reading records of file

printf("close files and exit \n");
hisclose_c(tno);
uvclose_c(tno_in);
uvclose_c(tno);
gettimeofday(&t_one,0);
elapsed=timedifference_msec(t_zero,t_one);
printf("\npprs executed in %f milliseconds.\n",elapsed);

}  // end of main
