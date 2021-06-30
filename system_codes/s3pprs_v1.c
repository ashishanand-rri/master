/*= S3PPRS - PreProcess the saras 3 data file */
/*& Ravi */
/*: Utility */
/*+
 S3PPRS* reads the miriad uv dataset from SARAS 3
 8193 channel data with 16 spectra per switch setting
 averages the 16 records per source name
 ignores pre-existing flags
 rejects RFI using Hampel filtering
 recomputes the LST
*/
/*@ vis */
/*         Name of the uvdataset; no defaults */
/*@ out */
/*         Name of the output uvdataset; no defaults */
/*@ tol */
/*         Values beyond tol*1.4826*MAD are rejected in the averaging; default 2.0 */
/*@ callst */
/*         Recompute lst if requested.  1=recompute, 0=do nothing (default) */
/*@ sitelatitude */
/*	       Site latitude (degrees); default 14.24232777778 */
/*@ sitelongitude */
/*		   Site longitude (degrees); default 77.612605555556 */
/*-- */

/*
c***************************************************************
c
c History
c
c Ravi 14Aug17 original version from s2pprs_v7.c
c
c***************************************************************
*/

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <maxdimc.h>
#include <miriad.h>
#include <sys/time.h>

#define PI 3.14159265
#define DATASIZE 8193

#define TOL 2.0
#define HW 12
#define BETA 61035.156250
#define INT_SPEC 0.067108864

void hampel(float darray[], int flags[], int num_chan, float tol, int hw);
double cal_lst(double utc_julian_date, double longitude_radians);




float timedifference_msec(struct timeval t0,struct timeval t1)
{
	return(t1.tv_sec-t0.tv_sec)*1000.0f + (t1.tv_usec-t0.tv_usec)/1000.0f;
}

int main(argc,argv)
int argc;
char *argv[];
{
	
struct timeval t0;
struct timeval t1;
float elapsed;
gettimeofday(&t0,0);

	
const char version[30]="s3pprs: version 1 14Aug17";
char visfile[50],outfile[50];
char status[30]="old";

int i,j,g,indx,k;
int tno,tno_in;
int ndata,nread,nread_store,nspect;
int datasize;
double preamble_store[3][4], preamble_auto11[4], preamble_auto22[4], preamble_cros12[4];
int hw,callst;
float tol;
int rflags[DATASIZE],iflags[DATASIZE], flags11[DATASIZE], flags22[DATASIZE];
float tarray[DATASIZE];
double gooddata,baddata;

int newstart,naverage[3];
int ant_code,ant1,ant2;
int newsource;
char oldsource[30];
char flag_name[30];
double time_var_store[1];
float data_array[3][16][2*DATASIZE];
int flags_array[16][DATASIZE];
float avdata[3][2*DATASIZE];
float avamp[2*DATASIZE];
float auto11[2*DATASIZE],auto22[2*DATASIZE],cros12[2*DATASIZE];

// Miriad-related variables
char source_name[30];
int flags[DATASIZE];
float flags_d[2*DATASIZE];
double preamble[4];
float data[2*DATASIZE];
char velocity_type[9] = {'V','E','L','O','-','O','B','S','\0'};
int ispec;
int fname_length;
int cnumber;
int var_ivalue[1];
static int npoints,nn;
int antenna1,antenna2,antenna3;
int access_result;
float baseline_value[1];
double p1,p2,p3;
double ant_az[3],ant_el[3];
double freq_channel1[1];
double freq_inc[1];
double time_var[1];
double coord_var[2];
double site_latitude[1],site_longitude[1];
float var_veldop[1];
float var_vsource[1];
double var_restfreq[1];
double sra[1],sdec[1];
double lo1[1],lo2[1],freq[1],freqif[1];
int mount[1];
float evector[1];
int nnpol,npol[1];
float jyperk[1],inttime[1],epoch[1];
double antpos[8];
float tpower[1];
float wfreq[1],wwidth[1];
int time_info[6][16];
int time_info1[6];
long int jjdd,yyyy,mmmm,dddd;
double julian_date;
double lst;
double longitude_radians,utc_julian_date;
float old_time,old_day,old_year;
float current_time,current_day,current_year;
double var, std_avg, std;
float sitelatitude, sitelongitude;
FILE *flag_out;

flag_name[0]='F';
flag_name[1]='L';
flag_name[2]='A';
flag_name[3]='G';
flag_name[4]='\0';

/* key calls to get user iniuts */
keyini_c(argc,argv);
keya_c("vis",visfile,"");
keya_c("out",outfile,"");
keyr_c("tol",&tol,TOL);
keyi_c("callst",&callst,0);
keyr_c("sitelatitude",&sitelatitude,13.44543);
keyr_c("sitelongitude",&sitelongitude,77.331982);
keyfin_c();

puts(version);

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
uvputvrd_c(tno,"antpos",antpos,8);
time_var[0]=18.0;
uvputvrd_c(tno,"time",time_var,1);
uvputvrd_c(tno,"lst",time_var,1);
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
site_latitude[0] = (double)(sitelatitude * PI/180.0);
uvputvrd_c(tno,"latitud",site_latitude,1);
site_longitude[0] = (double)(sitelongitude * PI/180.0);
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
hw = HW;
nn=30;
newstart=0;  // set to 1 when no longer first
gooddata=0.0;
baddata=0.0;

    for(;;)  // loop through all records in file 
    {
	uvread_c(tno_in,preamble,data,flags,ndata,&nread);
	if(nread == 0) break;

	/* preamble[3] is (double)(256*ant1+ant2) */
	ant_code=(int)preamble[3];
	ant2=ant_code%256;
	ant1=(ant_code-ant2)/256;
	
	if (ant1==1 && ant2==1) indx=0;
	else if (ant1==1 && ant2==2) indx=1;
	else if (ant1==2 && ant2==2) indx=2;

	// Get the source name and LST for this record
	uvgetvra_c(tno_in,"source",source_name,nn);
	uvgetvrd_c(tno_in,"lst",time_var,1);

	if (newstart == 0)  // if this is the first record ever
	{
	naverage[0]=0;
	naverage[1]=0;
	naverage[2]=0;
	
		for(i=0;i<ndata;++i) 
		{
			data_array[indx][naverage[indx]][2*i]=data[2*i];
			data_array[indx][naverage[indx]][2*i+1]=data[2*i+1];
		}
	time_var_store[0] = time_var[0];
	nread_store = nread;
	for(i=0;i<4;++i) preamble_store[indx][i]=preamble[i];
	for(i=0;i<30;++i) oldsource[i]=source_name[i];
	++naverage[indx];
	newstart = 1;
	}  // end of processing first record ever

	else  // got another record
	{								// if else processing non-first records
		newsource=0;
		for(i=0;i<30;++i)
		{
			if (source_name[i] != oldsource[i])
			 {
				newsource=1;
				break;
			}
		}
		if( newsource == 0)  // continuing on the oldsource
		{
			for(i=0;i<ndata;++i) 
			{
				data_array[indx][naverage[indx]][2*i]=data[2*i];
				data_array[indx][naverage[indx]][2*i+1]=data[2*i+1];
			}
			time_var_store[0] = time_var[0];
			nread_store = nread;
			for (i=0;i<4;++i) preamble_store[indx][i]=preamble[i];
			++naverage[indx];
		}
	
		if( newsource == 1)  // got a new source, implying old source records are all in place
		{
		// first process the older source
				
			if(naverage[0] != naverage[1] || naverage[1] != naverage[2] || naverage[2] != 16) 
			{ 
			printf(" Fault: num of recs for 11, 12 & 22 do not match! %d %d %d \n", 
				naverage[0],naverage[1],naverage[2]);
			}
		
			// 1. Hampel filter the 3*16 64 ms records
	
			for(j=0;j<naverage[indx];++j) 	// in the processing of the older source, 
											// naverage can have any indx value. 
			{
				for(i=0;i<ndata;++i) tarray[i]=data_array[1][j][2*i];  	// 12 real
				hampel(tarray, rflags, ndata, tol, hw);
				for(i=0;i<ndata;++i) tarray[i]=data_array[1][j][2*i+1];	// 12 imag
				hampel(tarray, iflags, ndata, tol, hw);
				for(i=0;i<ndata;++i) tarray[i]=data_array[0][j][2*i];	// 11 real
				hampel(tarray, flags11, ndata, tol, hw);
				for(i=0;i<ndata;++i) tarray[i]=data_array[2][j][2*i];	// 22 real
				hampel(tarray, flags22, ndata, tol, hw);
				
				for(i=0;i<ndata;++i) 
				{
					if (rflags[i]==0 || iflags[i]==0 || flags11[i]==0 || flags22[i]==0)
					{
						flags_array[j][i]=0;
					}
					else 
					{
						flags_array[j][i]=1;
					}
				}
			}
		
			// 2. Hampel filtering done, average the good data in the 16 64ms records 
	
			for(i=0;i<ndata;++i) // one channel at a time
			{
				nspect  = 0;    
				var     = 0.0;
				std_avg = 0.0;
				std     = 0.0; 
				avdata[0][2*i]=0.0;
				avdata[0][2*i+1]=0.0;
				avdata[1][2*i]=0.0;
				avdata[1][2*i+1]=0.0;
				avamp[2*i]=0.0;
				avamp[2*i+1]=0.0;
				avdata[2][2*i]=0.0;
				avdata[2][2*i+1]=0.0;
	
				for(j=0;j<naverage[indx];++j)  // loop over the 16 records
				{
					if (flags_array[j][i]==1)
					{
						avdata[0][2*i]+= data_array[0][j][2*i];
						avdata[0][2*i+1]+= data_array[0][j][2*i+1];
						avdata[1][2*i]+= data_array[1][j][2*i];
						avdata[1][2*i+1]+= data_array[1][j][2*i+1];
						avamp[2*i]+=sqrt(data_array[1][j][2*i]*data_array[1][j][2*i]+
										data_array[1][j][2*i+1]*data_array[1][j][2*i+1]);
						avdata[2][2*i]+= data_array[2][j][2*i];
						avdata[2][2*i+1]+= data_array[2][j][2*i+1];
						++nspect;
					}
				}
				if (nspect > 3)  // This scan of 16*64ms is considered good; so average the records
				{
					flags[i]=1;
					flags_d[2*i]=nspect;
					avdata[0][2*i] /= (float)nspect;
					avdata[0][2*i+1] /= (float)nspect;
					avdata[1][2*i] /= (float)nspect;
					avdata[1][2*i+1] /= (float)nspect;
					avamp[2*i] /= (float)nspect;
					avdata[2][2*i] /= (float)nspect;
					avdata[2][2*i+1] /= (float)nspect;
				}
				else 
				{
					flags[i]=0;
					flags_d[2*i]=0; 	// Flag the entire scan if number of good 64ms recs are < 4
					flags_d[2*i+1]=0;
				}
	
				if (flags[i]==1)	// If the scan is good, parse through all unflagged 64ms records
									// for the given channel and compute rms in 12imag
				{
					for (k=0; k<naverage[indx]; ++k)
					{	
						if (flags_array[k][i]==1)
//							var = var + (double)powf((data_array[1][k][2*i+1] - avdata[1][2*i+1]),2.0);
							var = var + (double)powf((sqrt(data_array[1][k][2*i]*data_array[1][k][2*i]+
										data_array[1][k][2*i+1]*data_array[1][k][2*i+1]) 
										- avamp[2*i]),2.0);

					}	
				var /= (float)(nspect-1);				
				std = sqrt(var);
				}
			flags_d[2*i+1]=std; // Enters computed std for unflagged points, else it is zero
			}  // end of one channel at a time
	
			// 3. Write the averaged spectra back into miriad	
	
			if (callst == 1)
			{
				utc_julian_date=preamble_store[indx][2]; 	// again, indx here can be anything 
										// from 0-2 as all LSTs are same
				longitude_radians = (double)site_longitude[0];
				lst=cal_lst(utc_julian_date, longitude_radians);
				time_var_store[0]=lst; 
			}
			uvputvr_c(tno,1,"source",oldsource,5);
			uvputvrd_c(tno,"lst",time_var_store,1);
		
			for (g=0; g < ndata; ++g)
			{	
				auto11[2*g]   = avdata[0][2*g];
				auto11[2*g+1] = avdata[0][2*g+1];
				auto22[2*g]   = avdata[2][2*g];
				auto22[2*g+1] = avdata[2][2*g+1];
				cros12[2*g]   = avdata[1][2*g];
				cros12[2*g+1] = avdata[1][2*g+1];
			}
			for(g=0; g<4; ++g) // this index of preamble has to be specific for auto11, auto22 and cross
			{
				preamble_auto11[g] = preamble_store[0][g];
				preamble_auto22[g] = preamble_store[2][g];
				preamble_cros12[g] = preamble_store[1][g];
			}
			uvwrite_c(tno,preamble_auto11,auto11,flags,nread_store);
			uvwrite_c(tno,preamble_cros12,cros12,flags,nread_store);
			uvwrite_c(tno,preamble_auto22,auto22,flags,nread_store);
		
			for(i=0; i<ndata; ++i)
			{
			if (flags[i]==1)
				gooddata+=flags_d[2*i];
			else
				baddata+=16;
			}
		
			uvputvr_c(tno,1,"source",flag_name,4);
			uvwrite_c(tno,preamble_cros12,flags_d,flags,nread_store);
		
			// Older source processed and written out
			// Start filling the arrays for the current source (first record)
		
			naverage[0]=0;//reset the averages
			naverage[1]=0;
			naverage[2]=0;
		
			for(i=0;i<ndata;++i) 
			{
				data_array[indx][naverage[indx]][2*i]=data[2*i];
				data_array[indx][naverage[indx]][2*i+1]=data[2*i+1];
			}
		
			time_var_store[0] = time_var[0];
			nread_store = nread;
			for (i=0;i<4;++i) preamble_store[indx][i]=preamble[i];
			for (i=0;i<30;++i) oldsource[i]=source_name[i];
			++naverage[indx];
		} // done with processing old plus new source
	
	} // end of if else processing a non-first record

	}  // end of loop through all records in file

	// flush the records already in buffers
 
		if(naverage[indx] > 0)
		{							// then process the last source of the file

			if(naverage[indx] != 16) 
			{ 
			printf(" Fault(at end): num of recs not 16! %d \n", naverage[indx]);
			}


			for(j=0;j<naverage[indx];++j)
			{
				for(i=0;i<ndata;++i) tarray[i]=data_array[1][j][2*i];
				hampel(tarray, rflags, ndata, tol, hw);
				for(i=0;i<ndata;++i) tarray[i]=data_array[1][j][2*i+1];
				hampel(tarray, iflags, ndata, tol, hw);
				for(i=0;i<ndata;++i) tarray[i]=data_array[0][j][2*i];
				hampel(tarray, flags11, ndata, tol, hw);
				for(i=0;i<ndata;++i) tarray[i]=data_array[2][j][2*i];
				hampel(tarray, flags22, ndata, tol, hw);

				for(i=0;i<ndata;++i) 
				{
					if (rflags[i]==0 || iflags[i]==0 || flags11[i]==0 || flags22[i]==0)
					{
						flags_array[j][i]=0;
					}
					else 
					{
						flags_array[j][i]=1;
					}
				}

			}		
	
		// Hampel filtering over, averaging starts
	
			for(i=0;i<ndata;++i)  // loop over channels
			{
			nspect  = 0;    
			var     = 0.0;
			std_avg = 0.0;
			std     = 0.0;
			avdata[0][2*i]=0.0;
			avdata[0][2*i+1]=0.0;
			avdata[1][2*i]=0.0;
			avdata[1][2*i+1]=0.0;
			avamp[2*i]=0.0;
			avamp[2*i+1]=0.0;
			avdata[2][2*i]=0.0;
			avdata[2][2*i+1]=0.0;

				for(j=0;j<naverage[indx];++j) 
				{
					if (flags_array[j][i]==1)
					{
					avdata[0][2*i]+= data_array[0][j][2*i];
					avdata[0][2*i+1]+= data_array[0][j][2*i+1];
					avdata[1][2*i]+= data_array[1][j][2*i];
					avdata[1][2*i+1]+= data_array[1][j][2*i+1];
					avamp[2*i]+=sqrt(data_array[1][j][2*i]*data_array[1][j][2*i]+
										data_array[1][j][2*i+1]*data_array[1][j][2*i+1]);
					avdata[2][2*i]+= data_array[2][j][2*i];
					avdata[2][2*i+1]+= data_array[2][j][2*i+1];
					++nspect;
					}
				}
				if (nspect > 3)
				{
				flags[i]=1;
				flags_d[2*i]=nspect;
				avdata[0][2*i] /= (float)nspect;
				avdata[0][2*i+1] /= (float)nspect;
				avdata[1][2*i] /= (float)nspect;
				avdata[1][2*i+1] /= (float)nspect;
				avamp[2*i] /= (float)nspect;
				avdata[2][2*i] /= (float)nspect;
				avdata[2][2*i+1] /= (float)nspect;
				}
				else 
				{
				flags_d[2*i]=0; // Flag the channel where number of averaged points are < 4
				flags_d[2*i+1]=0;
				flags[i]=0;
				}
				if (flags[i]==1)	// If the channel is good, parse through all
									// points for the given channel
									// to calculate the rms from unflagged points in 12imag
				{
					for (k=0; k<naverage[indx]; ++k)
					{	
						if (flags_array[k][i]==1)
//							var = var +  (double)powf((data_array[1][k][2*i+1] - avdata[1][2*i+1]),2.0); 	
							var = var + (double)powf((sqrt(data_array[1][k][2*i]*data_array[1][k][2*i]+
										data_array[1][k][2*i+1]*data_array[1][k][2*i+1]) 
										- avamp[2*i]),2.0);
					}	
				var /= (float)(nspect-1);				
				std = sqrt(var);
				}
			flags_d[2*i+1]=std; // Enters computed value for unflagged points, else it is zero
			}
		
			if (callst == 1)
			{
				utc_julian_date=preamble_store[indx][2];
				longitude_radians = (double)site_longitude[0];
				lst=cal_lst(utc_julian_date, longitude_radians);
				time_var_store[0]=lst; 
			}
		uvputvr_c(tno,1,"source",oldsource,5);
		uvputvrd_c(tno,"lst",time_var_store,1);
			for (g=0; g < ndata; ++g)
			{	
			auto11[2*g]   = avdata[0][2*g];
			auto11[2*g+1] = avdata[0][2*g+1];
			auto22[2*g]   = avdata[2][2*g];
			auto22[2*g+1] = avdata[2][2*g+1];
			cros12[2*g]   = avdata[1][2*g];
			cros12[2*g+1] = avdata[1][2*g+1];
			}
			for(g=0; g<4; ++g)
			{
			preamble_auto11[g] = preamble_store[0][g];
			preamble_auto22[g] = preamble_store[2][g];
			preamble_cros12[g] = preamble_store[1][g];
			}
		uvwrite_c(tno,preamble_auto11,auto11,flags,nread_store);
		uvwrite_c(tno,preamble_cros12,cros12,flags,nread_store);
		uvwrite_c(tno,preamble_auto22,auto22,flags,nread_store);
	
		for(i=0; i<ndata; ++i)
		{
		if (flags[i]==1)
			gooddata+=flags_d[2*i];
		else
			baddata+=16.0;
		}
	
		uvputvr_c(tno,1,"source",flag_name,4);
		uvwrite_c(tno,preamble_cros12,flags_d,flags,nread_store);
	
		} // end of processing the last source of the file

printf("Fraction of data bad %g\n",baddata/(gooddata+baddata));

printf("close files and exit\n");
hisclose_c(tno);
uvclose_c(tno_in);
uvclose_c(tno);
gettimeofday(&t1,0);
elapsed=timedifference_msec(t0,t1);
printf("\npprs executed in %f milliseconds.\n",elapsed);
} // end of main
 
