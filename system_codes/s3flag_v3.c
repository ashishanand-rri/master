/*= S3FLAG - Flagging the saras3 calibrated data file  */

/*& Ravi */
/*: Utility */
/*+
 S3FLAG reads the miriad uv dataset containing calibrated scans.
 Does a flagging of outliers.  Writes out a flagged file.
*/
/*@ vis */
/*         Name of the calibrated uvdataset; no defaults */
/*@ out */
/*         Name of the flagged uvdataset to write; no defaults */
/*@ chnl */
/*         start channel num (1-8193) to examine for RFI; default 1312 = 40 MHz */
/*@ chnh */
/*         end channel num (1-8193) to examine for RFI; default 2950 = 90 MHz */
/*@ chnm */
/*         mid channel num (2-8192) to use as pivot point for Taylor; default 2132 */
/*@ polyorder */
/*         Polynomial order to use for rebaseline step; default 10 */
/*@ htoll1 */
/*         Hampel tolerance value; values beyond tol*1.4826*MAD are rejected; default 5 */
/*@ htoll2 */
/*         Hampel tolerance value; values beyond tol*1.4826*MAD are rejected; default 3 */
/*@ maxstd */
/*         Maximum standard deviation acceptable ; default 6.0 */
/*@ autoh */
/*         Maximum acceptable value for auto correlation spectra; default 1500.0 */
/*@ autol */
/*         Minimum acceptable value for auto correlation spectra; default 100.0 */
/*@ crossh */
/*         Maximum acceptable value for cross correlation spectra; default 100.0 */
/*@ crossl */
/*         Minimum acceptable value for cross correlation spectra; default -300 */
/*-- */
/*

c***************************************************************
c
c History for s2flag
c
c Ravi 24Oct16 original version
c Saurabh 12Nov16 changed to include ACFs 11 and 22
c Saurabh 24Dec16 changed to include L1 and L2 flagging, some other minor changes
c Saurabh 26Dec16 changed to avoid fitting autocorrelations
c
c History for s3flag
c
c Ravi 27Aug17 original version from s2flag_v8
c Ravi/Saurabh  v2 is for low band 40-120 MHz
c Ravi 22Dec19  v3 is for ultra low band 40 - 90 MHz
c
c***************************************************************
*/

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <maxdimc.h>
#include <nrutil.h>
#include <miriad.h>

#define PI 3.14159265
#define DATASIZE 8193

#define FTOL 1.0e-3  /* fractional tolerance criterion for exiting fitting optimization */
#define EPS 1.0e-1   /* amoeba vertices are generated by incrementing coeffs with this value */

int NP, MP, NNPP, MMPP, chnm, indx;

int amoeba(float **p, float y[], int ndim, float ftol, float (*funk)(float []), int *nfunk);
float aelect_nr(unsigned long k, unsigned long n, float arr[]);
int npoints,ngoodpoints;
double total_points, flagged_points; 
float *xarray, *yarray, *yarray_pass, *yarray_12_res;
float *yarray_11, *yflags, *yflags_11,*yarray_12,*yflags_12,*yarray_22,*yflags_22 ;
int indx;
float select_nr(unsigned long k, unsigned long n, float arr[]);

/* polyval is called while determining fit y-values for subtraction from data */
float polyval(float xx,float ac[])
{
	int jj;
	float yy;

	xx = xx - (float)chnm;
	
	yy = 0.0;
	for(jj=1;jj<=NP;++jj) /* loop through all polynomial coefficients */
	{
	yy += ac[jj]*(powf(xx,(float)(jj-1)));
	}
	return (yy);
}		

/* funk is called by amoeba to get chisq for each guess coeff array */
float funk_L1(float ad[])
{
	int ii,jj;
	float xx,yy;
	float dev;

	dev=0.0;

	for(ii=1;ii<=npoints;++ii) /* loop through all spectral channels */
	{
		if (yflags[ii] == 1)
		{
		xx=xarray[ii] - (float)chnm;
		yy = 0.0;
			for(jj=1;jj<=NNPP;++jj) /* loop through all polynomial coefficients */
			{
			yy += ad[jj]*(powf(xx,(float)(jj-1)));
			}
		dev += fabs(yarray[ii] - yy);  // L1 minimization 
		}
	}
	return (dev);
}

float funk_L2(float ad[])
{
	int ii,jj;
	float xx,yy;
	float dev,dev2;

	dev2=0.0;

	for(ii=1;ii<=npoints;++ii) /* loop through all spectral channels */
	{
		if (yflags[ii] == 1)
		{
		xx=xarray[ii] - (float)chnm;
		yy = 0.0;
			for(jj=1;jj<=NNPP;++jj) /* loop through all polynomial coefficients */
			{
			yy += ad[jj]*(powf(xx,(float)(jj-1)));
			}
		dev = fabs(yarray[ii] - yy);
		dev2 += dev*dev;  // L2 minimization
		}
	}
	return (dev2);
}

/*
Used in funk_flag_zero
inspired from Python append function
*/
void append(float a[], int new_item, int N)
{
	int i;
	for (i=1;i<=N; i++)
	{
		if (a[i] == -1)
			break;	
	}
	a[i]=new_item;	
}	
	
/*
Used in funk_flag_zero
inspired from Python range function
*/
void range(float a[], int start, int end, int N)
{
	int i;
	for(i=start; i<end; i++)
		append(a, i, N);			
}

/*
Function used at the end of second round of L1 fitting
Takes in the data array and flag table
Modifies the flag table after flagging channels till two zero crossings are found
*/
void funk_flag_zero(float y[], float flag_arr[], int N)
{
	float *bad_pts_index;
	bad_pts_index = vector(1, N) ; //Array for storing the flagged channel indices of the current run	
	int num_zero_for, num_zero_bac, count;
	int i,j,k;

	for(k=1; k<=N; k++)
		bad_pts_index[k]=-1; //All points are initialized to -1

	for (i=1; i<=N; i++)
	{
		if(flag_arr[i]==0)
		{
			num_zero_for = 0;
			num_zero_bac = 0;
					
			//Backward run from the flagged channel
			for (j=i-1; j>=2; j--)
			{
				if(flag_arr[j]==0) // If any flagged channel is found on the way, skip
					continue;

				if (y[j]*y[j-1] < 0)  // If a zero crossing is found
				{
					num_zero_bac = num_zero_bac + 1; // Increment the counter
					if ((num_zero_bac == 2)||(j==2 && num_zero_bac < 2)) 
													// If zero crossing reaches two or
											     	// we reach beginning of spectra
					{
						// APPEND TO ARRAY
						range(bad_pts_index, j-1, i, N);  						
						num_zero_bac = 0;
						break;
					}
				}
			}

			//Forward
			for (j=i+1; j<=(N-1); j++)
			{
				if(flag_arr[j]==0)
					continue;
				
				if (y[j]*y[j+1] < 0) 
				{
					num_zero_for = num_zero_for + 1;
					if ((num_zero_for == 2)||(j==(N-2) && num_zero_for < 2))
					{
						//APPEND TO ARRAY
						range(bad_pts_index, i+1, j+2, N);  //If zero crossing reaches two or
															//we reach end of spectra				
						num_zero_for = 0;
						break;
					}
				}
			}
		}
	}
	for(k=1; k <= N; k++)
	{
		if (bad_pts_index[k]==-1) 
			break;
		else
			flag_arr[(int)bad_pts_index[k]]=0;			
	}	
}

/*
Computes standard deviation on the residuals (unflagged points)
*/
float std_dev(float y[], float f[], int N)
{
	int i, num1;
	float std1, sum1, mean1, var1;

	sum1=0;
	var1=0;
	num1=0;

	for(i=1; i<=N; i++)
	{
		if(f[i]==1)
		{
			sum1 = sum1 + y[i];
			num1 = num1 + 1;
		}
	}
	if (num1==0)
		return 100;
	
	mean1 = sum1/(float)num1;

	for(i=1; i<=N; i++)
	{
		if(f[i]==1)
			var1 = var1 + powf((y[i]-mean1),2);
	}
	var1 = var1/(float)num1;
	
	return sqrt(var1);
}

void clip(float spec_data[], float flag_data[], int N, int tag, float autoh1, 
											float autol1, float crossh1, float crossl1)
{
	int i;
	if(tag==0 || tag==2) // Autocorrelation spectra
	{
		for (i=1; i<=N; i++)
		{
			if (spec_data[i] > autoh1 || spec_data[i] < autol1)
				flag_data[i]=0;	
		}
	}
	if(tag==1)
	{
		for (i=1; i<=N; i++) // Cross correlation spectra
		{
			if (spec_data[i] > crossh1 || spec_data[i] < crossl1)
				flag_data[i]=0;	
		}
	}
}


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
    const char version[30]="s3flag: version 1 03Oct17";
    char visfile[50],outfile[50];
    char status[30]="old";
    int i,j,ii,nn, g, write_rec ;
    int tno,tno_in;
    int ndata,nread;
	int iwin,datasize;
    unsigned long ncal0,nobs0,ncal1,nobs1,nspect,ispect;
	int phasesw;
	 double preamble_store[3][4], preamble_11[4], preamble_22[4], preamble_12[4];
    float tcal, std, autoh, autol, crossh, crossl;
    float spect_array[3][2*DATASIZE];
    float auto_11[2*DATASIZE], auto_22[2*DATASIZE], cros_12[2*DATASIZE] ;
	float spect_real[3][DATASIZE],spect_real_randomized[DATASIZE] ;
	int spect_real_indices[DATASIZE];
    int flags_store[3][DATASIZE];
	double lst_store;

	// Miriad-related variables
    int flags[DATASIZE];
    float data[2*DATASIZE];
    double preamble[4];
	char source_name[30];
	char velocity_type[9] = {'V','E','L','O','-','O','B','S','\0'};
	int ispec;
	int fname_length;
	int cnumber;
	int var_ivalue[1];
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
	unsigned long pts;
	float tpower[1];
	float wfreq[1],wwidth[1];
	int time_info[6][16];
	int time_info1[6];
    long int jjdd,yyyy,mmmm,dddd;
    double julian_date;
	double lst;
	float old_time,old_day,old_year;
	float current_time,current_day,current_year;

	int ant1,ant2,ant_code;
	float maxstd;
	int chnl,chnh,polyorder,htol_l1, htol_l2;
	float thr[3];
	double datacount, flagcount;

	int num_rec;
    float *a_amoeba, **p_amoeba, *y_amoeba, ftol_amoeba, *fit_coeffs;
    int ndim_amoeba, nfunc_amoeba;
	float eps;
	float xx,yy;

    /* key calls to get user inputs */
    keyini_c(argc,argv);
    keya_c("vis",visfile,"");
    keya_c("out",outfile,"");
    keyi_c("chnl",&chnl,1312);
    keyi_c("chnh",&chnh,2950);
    keyi_c("chnm",&chnm,2132);
    keyi_c("polyorder",&polyorder,10);
    keyi_c("htoll1",&htol_l1,3);
    keyi_c("htoll2",&htol_l2,2);
    keyr_c("maxstd",&maxstd,18);
    keyr_c("autoh",&autoh,5000);
    keyr_c("autol",&autol,50);
    keyr_c("crossh",&crossh,200);
    keyr_c("crossl",&crossl,-200);
    keyfin_c();

	chnl -= 1;
	chnh -= 1;
	chnm -= 1;

    puts(version);

	if(chnl >= chnh) 
		{printf(" ....... min channel exceeds max channel: exiting \n"); return(0);}
	if(chnm < chnl || chnm > chnh)
		{printf(" ....... pivot channel is outside channel range: exiting \n"); return(0);}

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
	freq_inc[0]=  3.0517578125e-05; /* GHz */
	uvputvrd_c(tno,"sdf",freq_inc,1);
	site_latitude[0] = (double)(14.24232777778 * PI/180.0);
	uvputvrd_c(tno,"latitud",site_latitude,1);
	site_longitude[0] = (double)(77.612605555556 * PI/180.0);
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
	eps=EPS;

	NP = polyorder+1;
	MP = NP+1;
    a_amoeba = vector(1,NP);
    fit_coeffs = vector(1,NP);
    y_amoeba = vector(1,MP);
    p_amoeba = matrix(1,MP,1,NP);  
    ftol_amoeba=FTOL;

	xarray = vector(1,DATASIZE);
	yarray = vector(1,DATASIZE);
	yflags = vector(1,DATASIZE);
	yarray_11 = vector(1,DATASIZE);
	yflags_11 = vector(1,DATASIZE);
	yarray_22 = vector(1,DATASIZE);
	yflags_22 = vector(1,DATASIZE);
	yarray_12 = vector(1,DATASIZE);
	yflags_12 = vector(1,DATASIZE);
	yarray_pass = vector(1,DATASIZE);
	yarray_12_res = vector(1,DATASIZE);

	datacount=0.0;
	flagcount=0.0;
	total_points=0.0;
	flagged_points=0.0;
	num_rec=0;

    for(;;)
    {
	uvread_c(tno_in,preamble,data,flags,ndata,&nread);
	if(nread == 0) break;
	
	uvgetvra_c(tno_in,"source",source_name,nn);
	uvgetvrd_c(tno_in,"lst",time_var,1);

	if (source_name[2]=='C' && source_name[3]=='F') // No flagging for source named 'SPCF'
	{
		uvputvr_c(tno,1,"source",source_name,4);
		uvputvrd_c(tno,"lst",time_var,1);
		uvwrite_c(tno,preamble,data,flags,nread);
		continue; // Go to the next record
	}
	
	/* preamble[3] is (double)(256*ant1+ant2) */
	ant_code=(int)preamble[3];
	ant2=ant_code%256;
	ant1=(ant_code-ant2)/256;
	
	if (ant1==1 && ant2==1) {indx=0;num_rec+=1;}
	else if (ant1==1 && ant2==2) {indx=1;num_rec+=1;}
	else if (ant1==2 && ant2==2) {indx=2;num_rec+=1;}
	
	for(i=0;i<4;++i) preamble_store[indx][i]=preamble[i];

	lst_store = time_var[0];

	for(i=0;i<DATASIZE;++i) 
	{
	flags_store[indx][i]=flags[i];				/* Put flags into flags_store */
	}

	for(i=0;i<2*DATASIZE;++i)
	{
	spect_array[indx][i]=data[i];      			/* Put data into spect_array */
	}

	for(i=0;i<DATASIZE;++i) 
	{
	spect_real[indx][i]=spect_array[indx][2*i];    		/* Put real data into spect_real */
	}

	/* extract a segment */
	npoints = 0;
	for(i=chnl;i<chnh+1;++i)
	{
	++npoints;
	xarray[npoints]=(float)i;
	yarray[npoints]=spect_real[indx][i];     		/* Put segment of real data into yarray */
	yflags[npoints]=flags_store[indx][i];			/* Put corresponding segment of flags into yflags */
	}
	
    clip(yarray, yflags, npoints, indx, autoh, autol, crossh, crossl); 

/*
All 3 records, Auto 11, Auto 22 and Cross 12 are clipped. 
Hence the flags for auto 11/22 might not be exactly equal to cross 12
*/
	if (indx == 0)
	{
		for(i=1; i<=npoints; i++)
			{yflags_11[i] = yflags[i];			
			 yarray_11[i] = yarray[i];}	//data and flags stored in respective arrays		
	}
	if (indx == 1)
	{
		for(i=1; i<=npoints; i++)
			{yflags_12[i] = yflags[i];			
			 yarray_12[i] = yarray[i];}			
	}
	if (indx == 2)
	{
		for(i=1; i<=npoints; i++)
			{yflags_22[i] = yflags[i];			
			 yarray_22[i] = yarray[i];}			
	}
	
	if (indx == 1) //Proceed to this section only if it is cross correlation
	{
	/* set up amoeba and fit polynomial */

	ngoodpoints=0;
	for (i=1;i<=npoints;++i)
	{ 
		if(yflags[i]==1)
		{
		fit_coeffs[1] += yarray[i];
		++ngoodpoints;						/* ngoodpoints is number of unflagged pts in segment */ 
		}
	}
	datacount += (float)ngoodpoints;
	
	if(ngoodpoints > NP)					/* if sufficient good pts in segment, do a fit */
	{
		fit_coeffs[1] /= (float)ngoodpoints;  /* mean of the segment */
		
		/* Iterate to successively approximate a Taylor fit */
	
		for (ii=1; ii<=NP; ++ii)
		{
		NNPP = ii;
		MMPP = ii+1;
		ndim_amoeba=NNPP;
	
			for(i=1;i<=MMPP;++i)
			{
				for(j=1;j<=NNPP;++j)
				{
				p_amoeba[i][j]=fit_coeffs[j];
				if(i==(j+1)) p_amoeba[i][j] += eps+eps*fit_coeffs[j];
				}
			}
			for(i=1;i<=MMPP;++i)
			{
				for(j=1;j<=NNPP;++j)
				{
				a_amoeba[j]=p_amoeba[i][j];
				}
			y_amoeba[i]=funk_L1(a_amoeba);
			}
			
			amoeba(p_amoeba,y_amoeba,ndim_amoeba,ftol_amoeba,funk_L1,&nfunc_amoeba);
	
			for (i=1;i<=NNPP;++i)
			{
			fit_coeffs[i]=p_amoeba[1][i];		/* get the fit coefficients */
			}
		
		if (ii<NP) fit_coeffs[NNPP+1]=0.0;
	
		}
	
		for (i=1;i<=npoints;++i)
		{ 
		xx=xarray[i];
		yy=polyval(xx,fit_coeffs);
		yarray[i] -= yy;					/* subtract the fit from the segment */
		}
		pts = 0;			
		for (i=1;i<=npoints;++i)
		{
			if (yflags[i]==1)
			{
				pts++;	
				yarray_pass[pts] = fabs(yarray[i]);
			}
		}

		/* compute MAD for the segment */
		thr[indx] = select_nr(pts/2,pts,(yarray_pass));

		/* reject values exceeding htol * 1.4826 * MAD */
		thr[indx] *= 1.4826*htol_l1;

		for (i=1;i<=npoints;++i)
		{
			if(yflags[i]==1 && (fabs(yarray[i]) > fabs(thr[indx])))
			{
				 yflags[i]=0;
				 yflags_12[i]=0;
				 flagcount += 1.0;
			}
		}
	}

//L1, second round

	ngoodpoints=0;
	for (i=1;i<=npoints;++i)
	{ 
	yarray[i]=data[2*(chnl+i-1)];			/* reload the segment to refit updating flags*/
	if(yflags[i]==1)
		{
		fit_coeffs[1] += yarray[i];
		++ngoodpoints;						/* ngoodpoints is number of unflagged pts in segment */ 
		}
	}
	if(ngoodpoints > NP)					/* if sufficient good pts in segment, do a fit */
	{
		fit_coeffs[1] /= (float)ngoodpoints;  /* mean of the segment */
		
		/* Iterate to successively approximate a Taylor fit */
	
		for (ii=1; ii<=NP; ++ii)
		{
		NNPP = ii;
		MMPP = ii+1;
		ndim_amoeba=NNPP;
	
			for(i=1;i<=MMPP;++i)
			{
				for(j=1;j<=NNPP;++j)
				{
				p_amoeba[i][j]=fit_coeffs[j];
				if(i==(j+1)) p_amoeba[i][j] += eps+eps*fit_coeffs[j];
				}
			}
			for(i=1;i<=MMPP;++i)
			{
				for(j=1;j<=NNPP;++j)
				{
				a_amoeba[j]=p_amoeba[i][j];
				}
			y_amoeba[i]=funk_L1(a_amoeba);
			}
			
			amoeba(p_amoeba,y_amoeba,ndim_amoeba,ftol_amoeba,funk_L1,&nfunc_amoeba);
	
			for (i=1;i<=NNPP;++i)
			{
			fit_coeffs[i]=p_amoeba[1][i];		/* get the fit coefficients */
			}
		
		if (ii<NP) fit_coeffs[NNPP+1]=0.0;
	
		}
	
		for (i=1;i<=npoints;++i)
		{ 
		xx=xarray[i];
		yy=polyval(xx,fit_coeffs);
		yarray[i] -= yy;					/* subtract the fit from the segment */
		}
		pts = 0;			
		for (i=1;i<=npoints;++i)
		{
			if (yflags[i]==1)
			{
				pts++;	
				yarray_pass[pts] = fabs(yarray[i]);
			}
		}

		/* compute MAD for the segment */
		thr[indx] = select_nr(pts/2,pts,(yarray_pass));

		/* reject values exceeding htol * 1.4826 * MAD */
		thr[indx] *= 1.4826*htol_l1;

		for (i=1;i<=npoints;++i)
		{
			if(yflags[i]==1 && (fabs(yarray[i]) > fabs(thr[indx])))
				{
				 yflags[i]=0;
				 yflags_12[i]=0;
				 flagcount += 1.0;
				}
		}
		funk_flag_zero(yarray, yflags, npoints); // yarray contains the residuals
	}

//L2 fit
	/* And repeat the process to get a good RFI free fit */

	ngoodpoints=0;
	for (i=1;i<=npoints;++i)
	{ 
	yarray[i]=data[2*(chnl+i-1)];			/* reload the segment to fit again with updated flags */
	if(yflags[i]==1)
		{
		fit_coeffs[1] += yarray[i];
		++ngoodpoints;						/* ngoodpoints is number of unflagged pts in segment */ 
		}
	}

	if(ngoodpoints > NP)					/* if sufficient good pts in segment, do a fit */
	{
		fit_coeffs[1] /= (float)ngoodpoints;  /* mean of the segment */
		
		/* Iterate to successively approximate a Taylor fit */
	
		for (ii=1; ii<=NP; ++ii)
		{
		NNPP = ii;
		MMPP = ii+1;
		ndim_amoeba=NNPP;
	
			for(i=1;i<=MMPP;++i)
			{
				for(j=1;j<=NNPP;++j)
				{
				p_amoeba[i][j]=fit_coeffs[j];
				if(i==(j+1)) p_amoeba[i][j] += eps+eps*fit_coeffs[j];
				}
			}
			for(i=1;i<=MMPP;++i)
			{
				for(j=1;j<=NNPP;++j)
				{
				a_amoeba[j]=p_amoeba[i][j];
				}
			y_amoeba[i]=funk_L2(a_amoeba);
			}
			
			amoeba(p_amoeba,y_amoeba,ndim_amoeba,ftol_amoeba,funk_L2,&nfunc_amoeba);
	
			for (i=1;i<=NNPP;++i)
			{
			fit_coeffs[i]=p_amoeba[1][i];		/* get the fit coefficients */
			}
		
		if (ii<NP) fit_coeffs[NNPP+1]=0.0;
	
		}
		
		for (i=1;i<=npoints;++i)
		{ 
		xx=xarray[i];
		yy=polyval(xx,fit_coeffs);
		yarray[i] -= yy;					/* subtract the fit from the segment */
		}

		for (i=1; i<=npoints; ++i)
			yarray_12_res[i] = yarray[i];   /* create a residual array from this final L2 fit */

		pts = 0;			
		for (i=1;i<=npoints;++i)
		{
			if (yflags[i]==1)
			{
				pts++;	
				yarray_pass[pts] = fabs(yarray[i]);
			}
		}

		/* compute MAD for the segment */
		thr[indx] = select_nr(pts/2,pts,(yarray_pass));

		/* reject values exceeding htol * 1.4826 * MAD */
		thr[indx] *= 1.4826*htol_l2;
	
		for (i=1;i<=npoints;++i)
		{
			if(yflags[i]==1 && (fabs(yarray[i]) > fabs(thr[indx])))
			{
				 yflags[i]=0;
				 yflags_12[i]=0;
				 flagcount += 1.0;
			}
		}
	}
	}
	if (num_rec==3) // Process all the three recs and write all in the miriad uvfile
	{

		for (i=1;i<=npoints;++i)
		{
			if( (yflags_11[i]==0 || yflags_22[i]==0) && (yflags_12[i]==1) )
			{

			/*
			yflags_11, yflags_12 and yflags_22 changed once in this code during clipping.
			After that, only yflags array is changed for cross correlation spectra. 
			If a given channel is above threshold, we flag all the above three (+yflags). 
			Hence from only fitting perspective, all the above three should be same.
			However, since clipping is done independently for all 3 (and before fitting), 
			there might be some additional flagged channels in auto 11/22 which might not 
			be present in yflags_12
			
			It is to be noted that at a given time, cross is read first followed by auto 11/22. 
			Hence this equalization is performed 'after' fitting for cross spectra. 
			*/
				 yflags_12[i]=0;
				 flagcount += 1.0;
			}
		}	
	
		//condition of skipping the record
		std = std_dev(yarray_12_res, yflags_12, npoints);
		if (std > maxstd)  //maxstd being a parameter, default is 6.0 K
		{
			 printf(" Residual std: %g K compared to %g K \n",std,maxstd);
			 for (j = 1; j <=npoints ; ++j)
				yflags_12[j] = 0;
		}	
		
		for(i=chnl;i<=chnh;++i)
			flags[i] = yflags_12[i-chnl+1];			/* update flags */
		
		source_name[0]='S';
		source_name[1]='P';
		source_name[2]='C';
		source_name[3]='D';
		source_name[4]='\0';
		
		time_var[0] = lst_store;

		for(g=0; g < 2*ndata; ++g)
			{
			auto_11[g] = spect_array[0][g];
			auto_22[g] = spect_array[2][g];
			cros_12[g] = spect_array[1][g];
			}

		for(g=0; g < 4; ++g)
			{
			preamble_11[g] = preamble_store[0][g];
			preamble_22[g] = preamble_store[2][g];
			preamble_12[g] = preamble_store[1][g];
			}

		uvputvr_c(tno,1,"source",source_name,4);
		uvputvrd_c(tno,"lst",time_var,1);
		uvwrite_c(tno,preamble_12,cros_12,flags,nread);
		uvwrite_c(tno,preamble_11,auto_11,flags,nread);
		uvwrite_c(tno,preamble_22,auto_22,flags,nread);

		for(g=0; g<ndata; g++) //NOT npoints
			{
			if (flags[g]==0)
				flagged_points +=1;				
			}
		total_points += ndata;
		num_rec=0;
		
	}
	
    }

	free_vector (a_amoeba,1,NP);
	free_vector (fit_coeffs,1,NP);
	free_vector (y_amoeba,1,MP);
	free_matrix (p_amoeba,1,MP,1,NP);

	free_vector (xarray,1,DATASIZE);
	free_vector (yarray,1,DATASIZE);
	free_vector (yflags,1,DATASIZE);
	free_vector (yflags_11,1,DATASIZE);
	free_vector (yflags_12,1,DATASIZE);
	free_vector (yflags_22,1,DATASIZE);
	free_vector (yarray_11,1,DATASIZE);
	free_vector (yarray_12,1,DATASIZE);
	free_vector (yarray_22,1,DATASIZE);
	free_vector (yarray_pass,1,DATASIZE);

	printf("Num of chan pts examined %g, num flagged %g, fraction: %g \n", 
			datacount,flagcount,flagcount/datacount);
	printf("Total fraction of bad data :  %g \n", 
			flagged_points/total_points);

    printf("close files and exit\n");
    hisclose_c(tno);
    uvclose_c(tno_in);
    uvclose_c(tno);
    gettimeofday(&t_one,0);
    elapsed=timedifference_msec(t_zero,t_one);
    printf("\npprs executed in %f milliseconds.\n",elapsed);

}
