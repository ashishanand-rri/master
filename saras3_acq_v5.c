/* * saras3_acq_v4.c 8June 2018
 *
 * Acquire UTC from the GPS, 
 * and write out a MIRIAD dataset
 
 *	
 * 17Sep14 PC time use; GPS time not considered
 * 15Dec14 Setting switch states included
 * Change from v9 : removed sleep
 * 10June15 get time using only ftime; open/close output file for every successful data record
 *
 * 13Feb2018 : rename saras2_acq_vf2_16k.c to saras3_acq_v1.c
 *             change controls and state cycles from 4 to 6 to include
 *             "dicky" switch at front end, cal on/off, optical switch
 * 20Feb2018 : change from v1 to v2. v2 version is for three control lines.
 *             one for optical switch at base station.
 *             two for noise cal and dicky switch at antenna base.
 * 17Apr2018 : change from v2 to v3. 
		the new switch+LNA 

 * 8June2018 : change from v3 to v4. 
 * This program is used to acquire data from 2x8192 FFT firmware(saras2_v9_bnw_fft16k_2_Working6June2018_1648hrs.bin). 
 * FPGA firmware was modified from 8x2048 point FFt to 2x8192 point FFT. Earlier, the 8-point parallel FFT gave output 
 * in bit reversed order. In the 2x8192 point architecture, only one half of data(16384 points in the FFT) is brought 
 * out and is in normal order. In this code , at the place where data is stitched to obtain 16384 points, mapping is 
 * changed to reflect the normal order FFt output- 
 
*/
//
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <sys/ioctl.h>
#include <bits/socket.h>
#include <sys/select.h>
#include <netinet/in.h>
#include <netinet/ether.h>
#include <netinet/tcp.h>
#include <netinet/udp.h>
#include <netinet/ip.h>
#include <netinet/if_ether.h>
#include <net/ethernet.h>
#include <arpa/inet.h>
#include <netdb.h>
#include <signal.h>

#include <math.h>
#include <complex.h>

#include <time.h>
#include <sys/timeb.h>
#include <sys/time.h>
#include <maxdimc.h>
#include <miriad.h>

#include <termios.h>
#include <fcntl.h>
#include <sys/signal.h>
#include <complex.h> 


typedef char	        S8;  					/* Signed 8-bit integer (character) */
typedef unsigned char   U8;  					/* Unsigned 8-bit integer (byte)    */
typedef signed short    S16; 					/* Signed 16-bit integer (word)     */
typedef unsigned short  U16; 					/* Unsigned 16-bit integer (word)   */
typedef signed long     S32; 					/* Signed 32-bit integer (long)     */
typedef unsigned long   U32; 					/* Unsigned 32-bit integer (long)   */
typedef float           FLT; 					/* 4-byte single precision (float)  */
typedef double          DBL; 					/* 8-byte double precision (double) */

#define SIZE 256
#define PI 3.14159265358979

#define SOCKETS_IP_ADDRESS	"192.168.1.93"
#define SOCKETS_PORT		53503
#define SOCKETS_TIMEOUT		10

#define BAUDRATE B9600
#define device "/dev/mydev"

#define _POSIX_SOURCE 1
#define FALSE 0
#define TRUE 1
#define ESCAPE CTRL('[')
#define LOCAL_TIME_HR 5
#define LOCAL_TIME_MM 30

#define MSG_IN_ERROR     -1
#define MSG_IN_COMPLETE  0

#define TSIP_DLE         1
#define TSIP_IN_PARTIAL  2

#define MSG_IN_COMPLETE  0
#define TSIP_DLE         1
#define TSIP_IN_PARTIAL  2

#define DLE              0x10 					// TSIP packet start/end header         
#define ETX              0x03 					// TSIP data packet tail                
#define MAX_TSIP_PKT_LEN 300  					// max length of a TSIP packet 

void setup_gps();
void Parse0x8F (unsigned char ucData[], int nLen);

// System time related variables
const struct tm *tm ;
const struct stat st;
struct timeb cur_time_struct2;
struct tm *ut_time;
struct tm *cur_time_struct;
long cur_time_long;
unsigned short millitime;
float timesec[16];
time_t timep;
int time_usec;
clock_t startt,endt;

typedef unsigned long long u64;

u64 u64useconds;
struct timeval tv;



// UTC time variables    
int utc_hh,utc_mm,utc_ss,utc_day,utc_month,utc_year,utc_flag;

// Input-related variables
long int size,temp;
int a,b,c,d,e,x,y,z,count_spec=0;
const int M=8;                        				// Number of parallel, pipelined N-point FFTs 
const int N=2048;                     				// Length of each pipelined FFT
size_t count=0;
int tno,l,data_ind=1, flag=0, npkt=8192, dpkt, start=0, stop=0,len,next,next2;
long int zero_len=10, pktsize=1024, start_ind=0, stop_ind=0;

// Memory allocation variables
unsigned char *stage = NULL;
unsigned char *spectra_assuming_markers =NULL;
unsigned char *spectra_assuming_no_markers = NULL;

unsigned char *stage_new = NULL;
unsigned char *chan_data = NULL;
FILE *logfile, *dat_out;
unsigned int result;
int datasock; // Socket ids for various UDP logical links (sockets)

/* --- Function kbhit(), a keyboard lookahead monitor --- */
int kbhit(void)
{
  	int cnt = 0;
  	int error;
  	static struct termios Otty, Ntty;

  	tcgetattr(0, &Otty);
  	Ntty = Otty;

  	Ntty.c_iflag          = 0;       /* input mode                */
	Ntty.c_oflag          = 0;       /* output mode               */
  	Ntty.c_lflag         &= ~ICANON; /* raw mode */
  	Ntty.c_cc[VMIN]       = CMIN;    /* minimum characters to wait for  */
  	Ntty.c_cc[VTIME]      = CTIME;   /* minimum time to wait */

  	if (0 == (error = tcsetattr(0, TCSANOW, &Ntty))) 
	{
    		struct timeval      tv;
    		error     += ioctl(0, FIONREAD, &cnt);
    		error     += tcsetattr(0, TCSANOW, &Otty);
    		tv.tv_sec  = 0;
    		tv.tv_usec = 100;
    		select(1, NULL, NULL, NULL, &tv);
  	}

  	return (error == 0 ? cnt : -1 );
}

double cal_lst(double utc_julian_date, double longitude_radians)
{
	int julian_day,day_of_year;
	long int yyyy,mmmm,dddd,jjdd;
	double ut1_julian_date,time_ut1,gmst,lst,gast,dut1,utd,by,mjd;
	double dd,TT,HH,JD0,D0,DD;
	double ee, LL, OO, DP, eqeq;
	const double pi=3.14159265358979;

	mjd = utc_julian_date - 2400000.5;   /* Modified Julian Date */
	by = 1900.0 + (utc_julian_date - 2415020.31352) / 365.242198;   /* Besselian year */	

	utd = (0.022*sin(2*pi*by)) - (0.012*cos(2*pi*by)) - (0.006*sin(4*pi*by)) + (0.007*cos(4*pi*by)); 
	/* UTD = UT2-UT1 */ 	
	dut1 = 0.3523 - (0.00113*(mjd - 57955)) - utd;  /* DUT1 = UT1-UTC */

	ut1_julian_date = utc_julian_date + (double)dut1/(3600.0*24.0); 

	julian_day = (int)floor(ut1_julian_date + 0.5);
	time_ut1 = (ut1_julian_date+0.5) - (double)julian_day;
	if(time_ut1 >= 1.0) time_ut1 -= 1.0;
	if(time_ut1 < 0.0) time_ut1 += 1.0;
	JD0 = (double)julian_day - 0.5;
	HH = (ut1_julian_date - (double)JD0)*24.0;
	D0 = JD0 - 2451545.0;
	DD = ut1_julian_date - 2451545.0;
	TT = DD/36525.0;

	gmst = 6.697374558 + 0.06570982441908 * D0 + 1.00273790935 * HH + 0.000026 * TT*TT;
	while(gmst>=24.0) gmst -= 24.0;
	while(gmst<0.0) gmst += 24.0;

	ee = 23.4393 - 0.0000004 * DD;
	LL = 280.47 + 0.98565 * DD;
	OO = 125.04 - 0.052954 * DD;
	
	ee *= pi/180.0;
	LL *= pi/180.0;
	OO *= pi/180.0;

	DP = -0.000319*sin(OO) - 0.000024*sin(2*LL);
	eqeq = DP * cos(ee);

	gast = gmst + eqeq;

	lst = gast + (longitude_radians*180.0/pi)*(24.0/360.0);
	if(lst>24.0) lst -= 24.0;
	if(lst<0.0) lst += 24.0;

	lst *= (360.0*pi)/(24.0*180.0);
	
	return lst;
}

int main(int argc, char *argv[])
{ 
	setvbuf(stdout,NULL,_IONBF,0);
	int index;
	float  complex recieved_double[128];
	
	float  complex recieved_double_spectrum_11[16][16384];
	float  complex temp1[16384];
	
	float  complex recieved_double_spectrum_22[16][16384];
	float  complex temp2[16384];
	
	float  complex recieved_double_spectrum_12[16][16384];
	float  complex temp3[16384];
	float  complex tempcheck[16384];
	
	
	unsigned char *ptr1;
	unsigned char *ptr2;
	unsigned char *ptr3;
	unsigned char *ptr4;
	
	static float self11_float[16384];
	static float self22_float[16384];
	static float cross12_float[16384];
	
	
	float  complex  self11[16384];
	float  complex  self22[16384];
	float  complex  cross12[16384]; 
	
	//printf("checkpoint start of main");
  	long int i,k,j,ii,jj,zz,frame=0,n;
  	int SelfCorrPattern,CrossCorrPattern;
	int chan1_path[8] = {17,33,49,65,81,97,113,129};
	int chan2_path[8] = {18,34,50,66,82,98,114,130};
	int cc_real[8]    = {19,35,51,67,83,99,115,131};
	int cc_img[8]    = {20,36,52,68,84,100,116,132};
	int chan1_ind[8]  = {0,0,0,0,0,0,0,0};
	int chan2_ind[8] = {0,0,0,0,0,0,0,0};
	int cc_real_ind[8] = {0,0,0,0,0,0,0,0};
	int cc_img_ind[8] = {0,0,0,0,0,0,0,0};
	
	static long int Start_AutoCorr_Chan1[8][32774],Start_AutoCorr_Chan2[8][32774],
			Start_CrossCorrRealPath[8][32774],
			Start_CrossCorrImgPath[8][32774];
	static double selfcorr_chan1_complex[8][32774],selfcorr_chan2_complex[8][32774];
	static double real_crosscorr_complex[8][32774],img_crosscorr_complex[8][32774];
	static double before_reorder_selfcorr1[8][32774],reordered_selfcorr1[8][32774];
	static double before_reorder_selfcorr2[8][32774],reordered_selfcorr2[8][32774];
	static double before_reorder_real_crosscorr[8][32774],real_reordered_crosscorr[8][32774];
	static double before_reorder_img_crosscorr[8][32774],img_reordered_crosscorr[8][32774];
	static double final_selfcorr1[262152],final_selfcorr2[262152],final_real_crosscorr[262152],
			final_img_crosscorr[262152];
	static double selfcorr1_avg[262152],selfcorr2_avg[262152],crosscorr_ravg[262152],
			crosscorr_iavg[262152];
	double temp_sum,temp_selfcorr,real_crosscorr,img_crosscorr;

    	const char version[30]="saras3_acq_v1: 13FEB18";
    	char outfile[70],buf[20];
	char hostname[]="xport11";
	char command1[200],command2[200],command3[100];
	float temperature;
	temperature = 0;
	struct termios toptions;
	struct timeval timeout;      
	timeout.tv_sec = 5;
    	timeout.tv_usec = 0;
    	
	// Miriad-related variables
	char source_name[5];
	char velocity_type[9] = {'V','E','L','O','-','O','B','S','\0'};
	int ispec;
	int fname_length;
	int cnumber;
	int var_ivalue[1];
	static int npoints,nn;
	int antenna1,antenna2,antenna3;
	//static int flags[8193];
	static int flags[16384];
	
	int access_result;
	float baseline_value[1];
	double preamble[4];
	double p1,p2,p3;
	double ant_az[3],ant_el[3];
	double freq_channel1[1];
	double freq_inc[1];
	double time_var[1];
	double coord_var[2];
	double site_latitude[1],site_longitude[1],longitude_radians;
	static float data[16384];
	static float data2[32768];
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
	
	int data_size, header_size,sze=0,rcvbuf=16859136,oldbuf;//17907712
	int oldlen = sizeof(int);
    	struct sockaddr_in cliaddr,servaddr;
	socklen_t cliaddr_size;
    	int sock_raw;
	unsigned short iphdrlen;

	char ncal1_state,psw_state,optsw_state,pin;
	int inint,nint,nint4,fa,ct=1,p;
	char nint_buf[30];
	


	/* Arduino presets */
	//fa = open("/dev/ttyACM0",O_RDWR | O_NOCTTY);
	//if (fa == -1) 
	//{
		//perror("Unable to open port for Arduino\n");
		//return -1;
	//}

	//tcgetattr(fa, &toptions);
	//cfsetispeed(&toptions, B9600);
	//cfsetospeed(&toptions, B9600);
	//toptions.c_cflag &= ~PARENB;
	//toptions.c_cflag &= ~CSTOPB;
	//toptions.c_cflag &= ~CSIZE;
	//toptions.c_cflag |= CS8;
	//toptions.c_lflag |= ICANON;
	//tcsetattr(fa, TCSANOW, &toptions);

	/*query_gps();		//Query GPS for time.....commented
	time_info[0]=utc_hh;
	time_info[1]=utc_mm;
	time_info[2]=utc_ss;
	time_info[3]=utc_day;
	time_info[4]=utc_month;
	time_info[5]=utc_year;
	sprintf(outfile,"%04d-%02d-%02d_%02d%02d%02d.dat%c",
	time_info[5],time_info[4],time_info[3],
	time_info[0],time_info[1],time_info[2],'\0');*/
	
						//Query system time
	ftime(&cur_time_struct2);
	cur_time_long = cur_time_struct2.time;
  	cur_time_struct = (struct tm *) localtime(&cur_time_long);
	time_info1[0]=cur_time_struct->tm_hour;
	time_info1[1]=cur_time_struct->tm_min;
	time_info1[2]=cur_time_struct->tm_sec;
	time_info1[3]=cur_time_struct->tm_mday;
	time_info1[4]=cur_time_struct->tm_mon+1;
	time_info1[5]=cur_time_struct->tm_year + 1900;
	sprintf(outfile,"%04d-%02d-%02d_%02d%02d%02d.dat%c",
		  time_info1[5],time_info1[4],time_info1[3],
		  time_info1[0],time_info1[1],time_info1[2],'\0');
	//printf("Debug-checkpoint 1\n");	  
  		
	/* Setup parameters for writing the output miriad data file */

	/* set all flags to good */
	//for(i=0;i<8193;++i) flags[i]=1;
	for(i=0;i<16384;++i) flags[i]=1;

	npoints=8193; 
	
	
	
	//printf("Debug-checkpoint 2\n");
	/* Open output vis file */    	  	
	fprintf(stderr,"Open file with name: %s \n",outfile);
	
	//printf("Debug-checkpoint 3\n");
    
	uvopen_c(&tno,outfile,"new");
	hisopen_c(tno,"write");
	
	//printf("Debug-checkpoint 4\n");

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
	//var_ivalue[0]=8193;
	var_ivalue[0]=32768;
	uvputvri_c(tno,"nchan",var_ivalue,1);
	var_ivalue[0]=1;
	uvputvri_c(tno,"npol",var_ivalue,1);
	var_ivalue[0]=1;
	uvputvri_c(tno,"nspect",var_ivalue,1);
	var_ivalue[0]=-1;
	uvputvri_c(tno,"pol",var_ivalue,1);
	//var_ivalue[0]=8193;
	var_ivalue[0]=32768;
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
	//freq_inc[0]= 3.0517578125e-05; /* GHz */
	freq_inc[0]= 1;
	uvputvrd_c(tno,"sdf",freq_inc,1);
	site_latitude[0] = (double)(14.2423277777777778 * PI/180.0);
	uvputvrd_c(tno,"latitud",site_latitude,1);
	site_longitude[0] = (double)(77.61260555555555556 * PI/180.0);
	uvputvrd_c(tno,"longitu",site_longitude,1);
	antenna1=1;
	antenna2=2;
	baseline_value[0]=(float)(256*antenna1+antenna2);
	uvputvrr_c(tno,"baseline",baseline_value,1);
	coord_var[0]=0.0;
	coord_var[1]=0.0;
	uvputvrd_c(tno,"coord",coord_var,2);
	uvputvr_c(tno,1,"veltype",velocity_type,8);	

	// Close output files
	hisclose_c(tno);
	uvclose_c(tno);
	
	//printf("Debug-checkpoint 5\n");

	/* Allocate a buffer to hold the incoming UDP packet */
	unsigned char *buffer = (unsigned char *) malloc(1098);
	unsigned char *pkt_ctr = (unsigned char *) malloc(10);
	//printf("Debug-checkpoint 6\n");
	/* Allocate a buffer to hold the useful data from all UDP packets */
	if ((stage=calloc (npkt,pktsize)) == NULL)
  	{ perror ("calloc"); return -1; }
  	//printf("size of stage is %ld\n ", sizeof(spectra_assuming_no_markers));
  	
  	if ((stage_new=calloc (16380,768)) == NULL)
  	{ perror ("calloc"); return -1; }

  	
  	if ((spectra_assuming_markers=calloc (256,384)) == NULL)
  	{ perror ("calloc"); return -1; }
  	
  	
  	
  	if ((spectra_assuming_no_markers=calloc (256,1024)) == NULL)
  	{ perror ("calloc"); return -1; }
  	
  	//printf("size of spectra_assuming_no_markers below is %ld\n ", sizeof(spectra_assuming_no_markers));
  	
	printf("Debug-checkpoint 7\n");
	//* Allocate a buffer to hold headerless data */
	dpkt = npkt/2;
	if ((chan_data=calloc (dpkt,pktsize)) == NULL)
	{ perror ("calloc"); return -1; }
	//printf("Debug-checkpoint 8\n");
	strcpy(nint_buf,argv[1]);
	//printf("Debug-checkpoint 9\n");
  	sscanf(nint_buf,"%d",&nint4);
	//printf("Debug-checkpoint 10\n");
	fprintf(stderr,"Number of integrations (units of 6 states): %d \n",nint4);
	//printf("Debug-checkpoint 11\n");
	nint=nint4*6;
        
			// set initial state of CAL
			ncal1_state ='0';  
			pin = '4';
			write(fa,&pin,1); 
        		write(fa,&ncal1_state,1);

			// set initial state of dicky switch at antenna
			psw_state='0';
			pin = '5';
			write(fa,&pin,1); 
        		write(fa,&psw_state,1);

			// set initial state of opt switch at base station
			optsw_state='1';
			pin = '6';
			write(fa,&pin,1); 
        		write(fa,&optsw_state,1);
        		

	
	if ((datasock = socket(AF_PACKET, SOCK_RAW, htons(ETH_P_ALL))) < 0)
	  {
	   	perror("Socket Error");
	  
	  	return 1;
	  }
	
	setsockopt(datasock, SOL_SOCKET, SO_BINDTODEVICE, "eth0", strlen("eth0")+1);
	// Loop depending on the number of integrations
	for(inint=0;inint<nint;++inint)								
	{
	 //printf("Debug-checkpoint 12\n");
	  /* Open socket for communication and set options */



	   //printf("cliaddr_size is %d\n",cliaddr_size);
	   //printf("cliaddr_ip adress is %s\n",inet_ntoa(cliaddr.sin_addr));
	   //printf("cliaddr_port is is %hu\n",ntohs(cliaddr.sin_port));
	   //z=connect(datasock,(struct sockaddr *) &cliaddr,cliaddr_size);
	   //printf("Value of z is %d ",z);
	   //perror("connect");
	    //if (z=0)
	    //{
	      //printf("Socket connected\n");
	    //}
	    //else if(z<0)
	    //{
	      //printf("Socket NOT connected\n");
	    //}
	    
	    //else if(z>0)
	    //{
	      //printf("Socket NOT connected\n");
	    //}
	    
	    
	  //servaddr.sin_addr.s_addr=inet_addr("169.254.164.139");
	  //servaddr.sin_addr.s_addr=inet_addr("0.0.0.0");
	  //servaddr.sin_port=htons(5200);

	  
	  //a=setsockopt(datasock, SOL_SOCKET, SO_BINDTODEVICE, "eth0", strlen("eth0")+1);
	  //setsockopt(datasock, SOL_SOCKET, SO_RCVTIMEO, (char *)&timeout, sizeof(timeout));
	  //b=setsockopt(datasock, SOL_SOCKET, SO_RCVTIMEO, (char *)&timeout, sizeof(timeout));
	  //setsockopt(datasock, SOL_SOCKET, SO_RCVBUF, (char *)&rcvbuf, sizeof(int));
	  //c=setsockopt(datasock, SOL_SOCKET, SO_RCVBUF, (char *)&rcvbuf, sizeof(int));
	  
	  //printf("%d , %d , %d\n",a,b,c); 

          //if (inint%6==0)
	
	  /* Reading temperature from arduino.....currently commented*/
/* 	  tcflush(fa, TCIFLUSH);
	  n = read(fa, buf, 64);
	  buf[n] = 0;
	  temperature = atof(buf);
	 */

	// Set source name for this integration

				if(inint%6==0)
				{
					source_name[0]='O';
					source_name[1]='B';
					source_name[2]='S';
					source_name[3]='0';
					source_name[4]='0';
					source_name[5]='\0';
				}
				else if(inint%6==1)
				{
					source_name[0]='O';
					source_name[1]='B';
					source_name[2]='S';
					source_name[3]='1';
					source_name[4]='1';
					source_name[5]='\0';
				}	
				else if(inint%6==2)
				{
					source_name[0]='C';
					source_name[1]='A';
					source_name[2]='L';
					source_name[3]='0';
					source_name[4]='0';
					source_name[5]='\0';
				}	
				else if(inint%6==3)
				{
					source_name[0]='C';
					source_name[1]='A';
					source_name[2]='L';
					source_name[3]='0';
					source_name[4]='1';
					source_name[5]='\0';
				}	
				else if(inint%6==4)
				{
					source_name[0]='C';
					source_name[1]='A';
					source_name[2]='L';
					source_name[3]='1';
					source_name[4]='0';
					source_name[5]='\0';
				}	
				else if(inint%6==5)
				{
					source_name[0]='C';
					source_name[1]='A';
					source_name[2]='L';
					source_name[3]='1';
					source_name[4]='1';
					source_name[5]='\0';
				}	

	  //printf("Debug-checkpoint 14\n");
 
	  if (temperature < 40)
	  {
	  	if(inint%6==0)	// First integration of any cycle	
	    	{
			ncal1_state ='0';   
			pin = '4';
			write(fa,&pin,1); 
        		write(fa,&ncal1_state,1);

			psw_state='1';
			pin = '5';
			write(fa,&pin,1); 
        		write(fa,&psw_state,1);

			optsw_state='0';
			pin = '6';
			write(fa,&pin,1); 
        		write(fa,&optsw_state,1);
	    	}
	  	
	  	else if(inint%6 == 1) 
	  	{

			ncal1_state ='0'; 
			pin = '4';
			write(fa,&pin,1); 
        		write(fa,&ncal1_state,1);

			psw_state='1';
			pin = '5';
			write(fa,&pin,1); 
        		write(fa,&psw_state,1);

			optsw_state='1';
			pin = '6';
			write(fa,&pin,1); 
        		write(fa,&optsw_state,1);
	    	}

	  	else if(inint%6 == 2) 
	  	{

			ncal1_state ='0'; 
			pin = '4';
			write(fa,&pin,1); 
        		write(fa,&ncal1_state,1);

			psw_state='0';
			pin = '5';
			write(fa,&pin,1); 
        		write(fa,&psw_state,1);

			optsw_state='0';
			pin = '6';
			write(fa,&pin,1); 
        		write(fa,&optsw_state,1);
	    	}	

	  	else if(inint%6 == 3) 
	  	{

			ncal1_state ='0'; 
			pin = '4';
			write(fa,&pin,1); 
        		write(fa,&ncal1_state,1);

			psw_state='0';
			pin = '5';
			write(fa,&pin,1); 
        		write(fa,&psw_state,1);

			optsw_state='1';
			pin = '6';
			write(fa,&pin,1); 
        		write(fa,&optsw_state,1);
	    	}	

	  	else if(inint%6 == 4) 
	  	{

			ncal1_state ='1'; 
			pin = '4';
			write(fa,&pin,1); 
        		write(fa,&ncal1_state,1);

			psw_state='0';
			pin = '5';
			write(fa,&pin,1); 
        		write(fa,&psw_state,1);

			optsw_state='0';
			pin = '6';
			write(fa,&pin,1); 
        		write(fa,&optsw_state,1);
	    	}	

	  	else if(inint%6 == 5) 
	  	{

			ncal1_state ='1'; 
			pin = '4';
			write(fa,&pin,1); 
        		write(fa,&ncal1_state,1);

			psw_state='0';
			pin = '5';
			write(fa,&pin,1); 
        		write(fa,&psw_state,1);

			optsw_state='1';
			pin = '6';
			write(fa,&pin,1); 
        		write(fa,&optsw_state,1);
	    	}	

		printf(" %c%c%c%c%c Int %6d of %6d state:%2d %c %c %c\n",
			source_name[0],source_name[1],source_name[2],source_name[3],source_name[4],
			inint,nint,inint%6,ncal1_state,psw_state,optsw_state);
	 	//usleep(1000000);
		//usleep(100000);  // 100 ms sleep time to let switches settle
		//printf("Debug-checkpoint 15\n");
		/* Acquisition routine starts from here */
		memset (stage,0,sizeof(stage));
		memset (stage_new,0,sizeof(stage_new));
		memset (chan_data,0,sizeof(chan_data));
		//printf("Debug-checkpoint 16\n");
		// Start acquisition
		//startstop(2);
		//printf("Debug-checkpoint 17 in main accquisition\n");
		
		// ct is the packet counter which goes from 1 to 8192 per state; every 512 packets is one frame
		// We need to get 8192 packets of 1024 bytes
		// 16 frames of 4 spectra (11acf, 22acf, 12ccf_real, 12ccf_imag) of length 16384 channels
	       
		

	    //z=connect(datasock,(struct sockaddr *) &cliaddr,cliaddr_size);

	    while(ct <= 6144) 
		{
			//printf("1 here is segmentation fault\n");
			//printf("1");
			//printf("ct value is %d\n",ct);
			//memset (buffer,0,sizeof(buffer));
			//puts(buffer);			
			if (flag == 1)
			{
				flag = 0;
				memset (stage,0,sizeof(stage));
				memset (stage_new,0,sizeof(stage_new));
				memset (chan_data,0,sizeof(chan_data));
				//startstop(2);
				printf("Flag reset\n");
			}
		//printf("Debug-checkpoint 18\n");	
            // Receive a packet
			// Each recvfrom ought to give 1098 bytes, which is 1024 data bytes plus header
			if((data_size = recvfrom(datasock, buffer, 1066, 0,(struct sockaddr *) &cliaddr, 
				&cliaddr_size)) <= 0)
            {
            perror("\n Receive error");
	    //printf("data size is %d\n",data_size);
			//startstop(0);		
			ct = 1;
			flag = 1;
            }
            gettimeofday(&tv,NULL);
            u64useconds = (1000000*tv.tv_sec) + tv.tv_usec;
           //printf("2 here is segmentation fault\n");
            //char buffer[1098] = {0,0,0,0,0,1, };  //please uncomment just for test
	     //x=htons(cliaddr.sin_port);
	     //printf("===============================================\n");
	     //printf("data size is %d\n",data_size);
	     //printf("buffer size is %d\n",sizeof(buffer));
	     //printf("===============================================\n");
	     //printf("================BELOW IS DATA from FPGA SIDE\n");
	     //printf("%d",atoi(buffer));
	      //printf("%d",(int *)buffer);
	       //printf("%f\n",buffer);
	     //printf("\n===============================================\n");
         //printf("port %hu\n",x);
	     //printf("IP address %hu\n",cliaddr.sin_addr.s_addr);
	     //printf("===============================================\n");
	     
	     
			if(data_size == 1066)
			{	
			if((ct%512 == 0 || ct == 1 ))  // if at the beginning of a new frame
				{
				/*query_gps();		//Query GPS for time....currently commented
				time_info[0][frame]=utc_hh;
				time_info[1][frame]=utc_mm;
				time_info[2][frame]=utc_ss;
				time_info[3][frame]=utc_day;
				time_info[4][frame]=utc_month;
				time_info[5][frame]=utc_year;*/
									//Query system time.....currently used
				ftime(&cur_time_struct2);
				cur_time_long = cur_time_struct2.time;
				millitime=cur_time_struct2.millitm;
				cur_time_struct = (struct tm *) localtime(&cur_time_long);
				time_info[0][frame]=cur_time_struct->tm_hour;
				time_info[1][frame]=cur_time_struct->tm_min;
				time_info[2][frame]=cur_time_struct->tm_sec;
				time_info[3][frame]=cur_time_struct->tm_mday;
				time_info[4][frame]=cur_time_struct->tm_mon+1;
				time_info[5][frame]=cur_time_struct->tm_year + 1900;
				timesec[frame] = (float)cur_time_struct->tm_sec + (float)millitime/1000.0;

				frame++;
				}
				//printf("3 here is segmentation fault\n");
				next = (ct-1)*pktsize;  // N=1024
				next2= ct*pktsize;
			
				//Now process the current packet of 1024 bytes got in current recvfrom and store in buffer

				struct iphdr *iph = (struct iphdr*)(buffer + sizeof(struct ethhdr));
				
				struct udphdr *udph = (struct udphdr*)(buffer + iphdrlen  + sizeof(struct ethhdr));

                //printf("value of iphdr is %d",iph->protocol);
				switch (iph->protocol) //Check the Protocol and execute accordingly
				{
	                    
	                    //usleep(1000000);				
						case 17://UDP Protocol
						
							        //printf("Correct !!!! UDP packet encountered\n");
								//usleep(1000000);
								iphdrlen = iph->ihl*4;
								header_size =  sizeof(struct ethhdr) + iphdrlen + sizeof udph;// + 32;    // removing +32 as offset solved the problem of trailing 4 zeros
								//printf("\nHeader size is %d\n",header_size);
								memcpy (stage+next, buffer + header_size, pktsize);
								
								memcpy(recieved_double,buffer + header_size,pktsize);

								//printf("buffer+header_size",*(buffer+header_size));
								//for(int k=0;k<1024;k=k+8)   // This loop takes care of fact that out of 8 bytes last two bytes are not data ie markers 
								//{
								  //memcpy(stage_new ,stage+next+k,6);
							    //}
							    


								break;

						default://Other Protocols
								printf("Non-UDP packet encountered\n");

								//usleep(1000000);
								break;
				}
			 	ct = ct+1;
			 	
			 	//for (int k=0;k<128;k++)
			 	//{
				  //printf("recieved double is %f\n",recieved_double[k]);	
				//}
				//printf("ct value inside if is %d\n",ct);

				//usleep(1000000);
			}
			
			//printf("ct value outside any if: %d \n",ct);
			if(ct>128 &&(ct%128==1) && ((ct/128)%3)==1)    // 512 packets correspond to a frame when packet counter reaches 128 it must be end of self1 data , given 4 spectras are in 1 frame 
			{
				index=(ct-128)/(128*3);
				//printf("4 here is segmentation fault\n");
				//printf("Index self11 is %d\n",index);
                //printf("next2 value is %ld\n",next2);
                //printf("next value is %ld\n",next);
			    //printf("ct in self11 else : %d \n",ct);
				//memcpy(spectra_assuming_markers,stage_new,(256*768)/2);
				//printf("error possible here1");
				//memcpy(recieved_double_spectrum_11,stage+next-(128*pktsize),128*1024);
				//memcpy(&recieved_double_spectrum_11[index][16384],stage+next-(128*pktsize),128*1024);
				memcpy(temp1,stage+next-(127*pktsize),128*1024);
				//memcpy(tempcheck,stage,131072);
				//printf("error possible here2");
				for(int x=0;x<16384;x++)
				{
				  //printf("temp1 real,imag is %f+j%f\n",creal(temp1[x]),cimag(temp1[x]));
				 //printf("temp1 real,imag is %f+j%f\n",creal(tempcheck[x]),cimag(tempcheck[x]));
				 recieved_double_spectrum_11[index][x] =temp1[x];	
				}
			}
			
			
	         else if(ct>128 &&(ct%128==1) && ((ct/128)%3)==2)    // 512 packets correspond to a frame when packet counter reaches 128 it must be end of self1 data , given 4 spectras are in 1 frame 
			{
				//printf("5 here is segmentation fault\n");
				index=(ct-256)/(128*3);
				//printf("Index self22 is %d\n",index);
				//printf("next2 value is %ld\n",next2);
				//ptr2=stage+next;
			     //printf("ct in self22 else : %d \n",ct);
				//memcpy(spectra_assuming_markers,stage_new,(256*768)/2);
				//printf("error possible here1");
				//memcpy(&recieved_double_spectrum_22[index][16384],stage+next-(128*pktsize),128*1024);
				memcpy(temp2,stage+next-(127*pktsize),128*1024);
			        for(int y=0;y<16384;y++)
				{
				 recieved_double_spectrum_22[index][y] =temp2[y];	
				 //printf("temp2 real,imag is %f+j%f\n",creal(temp2[y]),cimag(temp2[y]));
				}
				//printf("error possible here2");
			}
			
		        else if(ct>128 &&(ct%128==1) && ((ct/128)%3)==0)     // 512 packets correspond to a frame when packet counter reaches 128 it must be end of self1 data , given 4 spectras are in 1 frame 
			{
				//printf("6 here is segmentation fault\n");
				index=(ct-384)/(128*3);
				//printf("Index cross12a is %d\n",index);
				//printf("next2 value is %ld\n",next2);
				//ptr3=stage+next;
			   // printf("ct in cross12a else : %d \n",ct);
				//memcpy(spectra_assuming_markers,stage_new,(256*768)/2);
				//printf("error possible here1");
				//memcpy(&recieved_double_spectrum_12[index][16384],stage+next-(128*pktsize),128*1024);
				
				memcpy(temp3,stage+next-(127*pktsize),128*1024);
			    for(int z=0;z<16384;z++)
				{
				  //printf("temp3 real,imag is %f+j%f\n",creal(temp3[z]),cimag(temp3[z]));	
				 recieved_double_spectrum_12[index][z] =temp3[z];	
				}
				
				
				//printf("error possible here2");
			}
			
			

			
		//printf("ct value inside if is %d\n",ct);
		}
		
		
		
		
		//printf("\n==============================================================\n");
		//for(int m=0;m<16;m++)
		//{
		//for (int k=0;k<16384;k++)
	       //{
		    ////printf(" index:[%d ,%d],self11value: %f+j%f\n",m,k,creal(recieved_double_spectrum_11[m][k]),cimag(recieved_double_spectrum_11[m][k]));	
		   //}
		
	    //}
		//printf("\n==============================================================\n");
		//for(int m=0;m<16;m++)
		//{
		//for (int k=0;k<16384;k++)
	       //{
		    ////printf(" index:[%d ,%d],self11value: %f+j%f\n",m,k,creal(recieved_double_spectrum_22[m][k]),cimag(recieved_double_spectrum_22[m][k]));	
		   //}
		
	    //}
		//printf("\n==============================================================\n");
		//for(int m=0;m<16;m++)
		//{
		//for (int k=0;k<16384;k++)
	       //{
		    ////printf(" index:[%d ,%d],self11value: %f+j%f\n",m,k,creal(recieved_double_spectrum_12[m][k]),cimag(recieved_double_spectrum_11[m][k]));	
		   //}
		
	    //}
	       //printf("\n==============================================================\n");
		

		//printf("Now outside of loop of ct<16380");
		//printf("size of spectra_assuming_no_markers below is %ld ", sizeof(spectra_assuming_no_markers));
		
    	// Got 16 frames for this switch setting. Stop acquisition.
		//startstop(0);
		ct = 1;
		
		
	


                if(kbhit() == 0)
                {

				// Open output files in append mode
				uvopen_c(&tno,outfile,"append");
				hisopen_c(tno,"append");

				// write the header record here with source name
				uvputvr_c(tno,1,"source",source_name,5);


                for(frame=0;frame<16;frame++)
                {		
                     //FILE *file1;
                     //FILE *file2;
                     //FILE *file3;
                     //printf("7 here is segmentation fault\n");
						//char buffer1[32]; // The filename buffer.
                    	
						//// Put "file" then k then ".txt" in to filename.
						//snprintf(buffer1, sizeof(char)*100, "/home/pi/Desktop/acquisistion_logs/frame%i_%i_self11.txt", frame,inint);
						//file1 = fopen(buffer1, "w");                     
                     
                    for(int x=0;x<16384;x++)
                    {
						self11[x]=recieved_double_spectrum_11[frame][x];

                        //fprintf(file1,"%f+j%f\n",creal(self11[x]),cimag(self11[x]));
							

						//printf(" (%d) %f+j%f ",x,creal(recieved_double_spectrum_11[frame][x]),cimag(recieved_double_spectrum_11[frame][x]));


                    }	
                    //fclose(file1);
                    
				
                   //printf("\n==============================================================\n");
                   						
                   //char buffer2[32]; // The filename buffer
                   //snprintf(buffer2, sizeof(char)*100, "/home/pi/Desktop/acquisistion_logs/frame%i_%i_cross12.txt", frame,inint);
                   //file2 = fopen(buffer2, "w");
                     for(int y=0;y<16384;y++)
                    {
						
						cross12[y]=recieved_double_spectrum_12[frame][y];
                        //fprintf(file2,"%f+j%f\n",creal(cross12[y]),cimag(cross12[y]));
						//printf("(%d) %f+j%f",y,creal(recieved_double_spectrum_12[frame][y]),cimag(recieved_double_spectrum_12[frame][y]));
						


                    }	
                    //fclose(file2);	
                    //printf("\n==============================================================\n");
                    

						
                    	//char buffer3[32]; // The filename buffer
						//snprintf(buffer3, sizeof(char)*100, "/home/pi/Desktop/acquisistion_logs/frame%i_%i_self22.txt", frame,inint);
						//file3 = fopen(buffer3, "w");				        
					
                    
                     for(int z=0;z<16384;z++)
                    {

						self22[z]=recieved_double_spectrum_22[frame][z];
						                    

                        //fprintf(file3,"%f+j%f\n",creal(self22[z]),cimag(self22[z]));	
						
						
						
						//printf("(%d) %f+j%f",z,creal(recieved_double_spectrum_22[frame][z]),cimag(recieved_double_spectrum_22[frame][z]));


                    }						
                    //fclose(file3);	
                   //printf("\n==============================================================\n");
                   

					preamble[0]=0.0;  /* u coordinate */
					preamble[1]=0.0;  /* v coordinate */
	
					// assemble the current JD in preamble[2] from the array time_info that contains UTC
					yyyy=time_info[5][frame];
					mmmm=time_info[4][frame];
					dddd=time_info[3][frame];
					jjdd = ( 1461 * ( yyyy + 4800 + ( mmmm - 14 ) / 12 ) ) / 4 +
							( 367 * ( mmmm - 2 - 12 * ( ( mmmm - 14 ) / 12 ) ) ) / 12 -
							( 3 * ( ( yyyy + 4900 + ( mmmm - 14 ) / 12 ) / 100 ) ) / 4 +
					dddd - 32075;
					julian_date = (double)jjdd -0.5;
					julian_date += (double)timesec[frame]/(3600.0*24.0) + 
							(double)time_info[1][frame]/(60.0*24.0) + (double)time_info[0][frame]/24.0;
					preamble[2] = julian_date;
	
					longitude_radians = (double)site_longitude[0];
									lst = cal_lst(julian_date, longitude_radians);
									time_var[0] = lst;
					
					
					//printf("Debug checkpoint before miriad write");				
					antenna1=1;
					antenna2=1;
					preamble[3]=(double)(256*antenna1+antenna2);

					// write the data record here
					//for(i=0;i<8193;++i) //Mayuri 27 Jun 2017
					//{
							////data[2*i]=(float)selfcorr1_avg[i];
							//data[2*i]=(float)chan_data[i];
							//data[2*i+1]=0.0;
					//}
					//printf("Debug checkpoint before miriad self2 write");
					
					//for(i=0;i<32768;++i) //Mayuri 27 Jun 2017
					//{
                     //data2[i]=spectra_assuming_no_markers[i];      
					//}
					
					
	
					uvputvrd_c(tno,"lst",time_var,1);
					//uvwrite_c(tno,preamble,data,flags,npoints);
					uvwrite_c(tno,preamble,self11,flags,16384);
	                 
					antenna1=2;	
					antenna2=2;
					preamble[3]=(double)(256*antenna1+antenna2);

					uvputvrd_c(tno,"lst",time_var,1);                                
					//uvwrite_c(tno,preamble,data,flags,npoints);
					uvwrite_c(tno,preamble,self22,flags,16384);
	                
					antenna1=1;
					antenna2=2;
					preamble[3]=(double)(256*antenna1+antenna2);

				
					uvputvrd_c(tno,"lst",time_var,1);
					uvwrite_c(tno,preamble,cross12,flags,16384);
					//uvwrite_c(tno,preamble,data,flags,npoints);

                } // End of all frames

				// Close output files at end of writing all 16 frame data
				hisclose_c(tno);
				uvclose_c(tno);

                } // End of kbhit == 0
		
				frame=0;
				memset(chan1_ind, 0, sizeof(chan1_ind));
				memset(chan2_ind, 0, sizeof(chan2_ind));
				memset(cc_real_ind, 0, sizeof(cc_real_ind));
				memset(cc_img_ind, 0, sizeof(cc_img_ind));
			
				if (kbhit() != 0) /* and if a keyboard hit is detected: exit from loop */
				{
				printf("\n Detected a keyboard hit..stopping..\n");

				// reset switches

			ncal1_state ='0'; 
			pin = '4';
			write(fa,&pin,1); 
        		write(fa,&ncal1_state,1);

			psw_state='0';
			pin = '5';
			write(fa,&pin,1); 
        		write(fa,&psw_state,1);

			optsw_state='0';
			pin = '6';
			write(fa,&pin,1); 
        		write(fa,&optsw_state,1);

	      
			//usleep(100000);  // 100 ms sleep time to let switches settle

				printf("  Counter reset \n");

				system(command3);
				free(buffer);
				free(stage);
				free(chan_data);
				close(fa);
				
				printf("  Miriad data written\n");
			
				exit(0);
				}
	  }  // End of code executed if arduino temp is ok....condition check commented at beginning of this loop

	// if the temperature measured by arduino is too high then close and quit....currently commented
	  else  
	  {

		// Switch reset and exit 

			ncal1_state ='0'; 
			pin = '4';
			write(fa,&pin,1); 
        		write(fa,&ncal1_state,1);

			psw_state='0';
			pin = '5';
			write(fa,&pin,1); 
        		write(fa,&psw_state,1);

			optsw_state='0';
			pin = '6';
			write(fa,&pin,1); 
        		write(fa,&optsw_state,1);
	      
			//usleep(100000);  // 100 ms sleep time to let switches settle

		system(command3);
		free(buffer);
		free(stage);
		free(chan_data);
		close(fa);
		
		printf("Miriad data written..counter reset..quitting because of high temp %f \n", temperature);
		exit(0);
	  }

	} // End of all integrations
	
	close(datasock);

			ncal1_state ='0';   // reset switches
			pin = '4';
			write(fa,&pin,1); 
        		write(fa,&ncal1_state,1);

			psw_state='0';
			pin = '5';
			write(fa,&pin,1); 
        		write(fa,&psw_state,1);

			optsw_state='0';
			pin = '6';
			write(fa,&pin,1); 
        		write(fa,&optsw_state,1);

			//usleep(100000);  // 100 ms sleep time to let counter settle

	system(command3);        
	free(buffer);
	free(stage);
	free(chan_data);
	close(fa);
	printf("Miriad data written,End of Code , Bye!!!\n");

	return 0;
}

