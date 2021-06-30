/*
*  hampel.c
*
*  Hampel filter spectra with threshold set at tol*1.4826*MAD
*
*  darray is the input spectrum
*  flags is the integer array of flags: flag is set to 0 if bad channel
*  num_chan is the number of channels in the array 
*  hw specifies the window size: window = 2*hw+1.
*  tol is the tolerance: values beyond tol*MAD are marked flag = 0
*  
*/

#include <math.h>
#include <nrutil.h>

#define DATASIZE 4097

float select_median(unsigned long k, unsigned long n, float arr[]);

int hampel(float darray[], int flags[], int num_chan, float tol, int hw)
{
int i, j;
int start_channel,end_channel;
int ihamp,jhamp;
int iratio;
float ratio1,ratio2,sk;
float ratio_array1[DATASIZE];
start_channel=hw;
end_channel= (num_chan-1)-hw;

// Reset all flags to good
for(i=0;i<num_chan;++i) flags[i]=1;

/* hampel filter of window size 2*hw+1 points */
	for(i=start_channel;i<=end_channel;++i)
	{
		ihamp=i-hw;
		jhamp=i+hw;
		iratio=0;
		for(j=ihamp;j<=jhamp;++j)
			{
			iratio += 1;
			ratio_array1[iratio]=darray[j];
			}
		ratio1=select_median((iratio+1)/2,iratio,ratio_array1);  // median
		iratio=0;
		for(j=ihamp;j<=jhamp;++j)
			{
			iratio += 1;
			ratio_array1[iratio]=fabs(darray[j]-ratio1);  // absolute deviations from median
			}
		ratio2=select_median((iratio+1)/2,iratio,ratio_array1);  // median absolute deviation

		sk=1.4826*ratio2;

		if( fabs(darray[i] - ratio1) > tol*sk ) flags[i]=0;

	}
return 0;
}

#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
float select_median(unsigned long k, unsigned long n, float arr[])
{
	unsigned long i,ir,j,l,mid;
	float a,temp;

	l=1;
	ir=n;
	for (;;) {
		if (ir <= l+1) {
			if (ir == l+1 && arr[ir] < arr[l]) {
				SWAP(arr[l],arr[ir])
			}
			return arr[k];
		} else {
			mid=(l+ir) >> 1;
			SWAP(arr[mid],arr[l+1])
			if (arr[l+1] > arr[ir]) {
				SWAP(arr[l+1],arr[ir])
			}
			if (arr[l] > arr[ir]) {
				SWAP(arr[l],arr[ir])
			}
			if (arr[l+1] > arr[l]) {
				SWAP(arr[l+1],arr[l])
			}
			i=l+1;
			j=ir;
			a=arr[l];
			for (;;) {
				do i++; while (arr[i] < a);
				do j--; while (arr[j] > a);
				if (j < i) break;
				SWAP(arr[i],arr[j])
			}
			arr[l]=arr[j];
			arr[j]=a;
			if (j >= k) ir=j-1;
			if (j <= k) l=i;
		}
	}
}
#undef SWAP


