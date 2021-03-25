
// Client side implementation of UDP client-server model 
#include <stdio.h> 
#include <stdlib.h> 
#include <unistd.h> 
#include <string.h> 
#include <sys/types.h> 
#include <sys/socket.h> 
#include <arpa/inet.h> 
#include <netinet/in.h> 
#include <time.h>
#include <unistd.h> 
#define PORT     4888 
#define MAXLINE 1024 
#include<complex.h>
#include<math.h>
  
// Driver code 
int main() 
{ 
    int sockfd;
    int z,counter=0;
    int j=csqrt(-1);

    char buffer[MAXLINE]; 
    float complex  array11[134*134];
    float complex  array22[134*134];
    float complex  array12[134*134];
    float complex  arraysend[134];
    
    FILE *file_ptr;

    for(int i=0;i<134*134;i++)
    {
     
     //if(i>3)
     //{
     //arrayOne[i]=i;
     //}
      //arrayOne[i].real=1;
     // arrayOne[i].imag=3;
     //creal(arrayOne[i]):=1;
     //cimag(arrayOne[i]):=2;
     array11[i]=CMPLX(1,0);   
    array12[i]=CMPLX(3,0);
    array22[i]=CMPLX(2,0);
     
     //array11[i]=CMPLX(i,0);   
     //array12[i]=CMPLX(i+1,0);
     //array22[i]=CMPLX(i+3,0);
     

     //printf("Array one element real and imag is %f,%f\n",creal(array11[i]),cimag(array11[i]));
    }   

    
    printf("size of one element of array 11 is %ld Bytes\n",sizeof array11);
    printf("size of one element of array 12 is %ld Bytes\n",sizeof array12);
    printf("size of one element of array 22 is %ld Bytes\n",sizeof array22);
    sleep(1);

    //printf("size of whole array is %ld\n",sizeof array11);
    struct sockaddr_in     servaddr; 
  
    // Creating socket file descriptor 
    if ( (sockfd = socket(AF_INET, SOCK_DGRAM, 0)) < 0 )
     { 
        perror("socket creation failed"); 
        exit(EXIT_FAILURE); 
     } 
   //if(sockfd=socket(AF_PACKET,SOCK_RAW,htons(ETH_P_ALL))<0)
   //{
     //perror("socket")
   // }  

    //memset(&servaddr, 0, sizeof(servaddr)); 
      
    // Filling server information 
    servaddr.sin_family = AF_INET; 
    servaddr.sin_port = htons(PORT); 
    servaddr.sin_addr.s_addr = inet_addr("169.254.164.139"); 
	//servaddr.sin_addr.s_addr = inet_addr("169.254.33.163");
          //setsockopt(sockfd, SOL_SOCKET, SO_BINDTODEVICE, "eno1", strlen("eth0")+1);
	  //a=setsockopt(datasock, SOL_SOCKET, SO_BINDTODEVICE, "eth0", strlen("eth0")+1);
	  //setsockopt(datasock, SOL_SOCKET, SO_RCVTIMEO, (char *)&timeout, sizeof(timeout));
	  //b=setsockopt(datasock, SOL_SOCKET, SO_RCVTIMEO, (char *)&timeout, sizeof(timeout));
	  //setsockopt(datasock, SOL_SOCKET, SO_RCVBUF, (char *)&rcvbuf, sizeof(int));
	  //c=setsockopt(datasock, SOL_SOCKET, SO_RCVBUF, (char *)&rcvbuf, sizeof(int));
         
    int n, len,frame=0;
   printf("Entering sending loop"); 
   while(frame<96)
{   
	
	if(counter==384)
	{
		counter=0;
		frame=frame+1;
		printf("\n==============================================================\n");
		printf("frame number is %d\n",frame);
		printf("\n==============================================================\n");
		
	}
	
	if ((frame+1)%16==0)
	{
		//usleep(182400);
		usleep(248400);
	}
     printf("Counter is %d\n",counter);
     //printf("\nsent packet till now is %d\n",counter);
     //sleep(1);
   // For self11 data
     if( floor(((counter/128)%3)==0))
    {
     for (int m=0;m<128;m=m+1)
     {
        //printf("\nm+128*counter %d , counter is %d , m is %d\n",m+128*counter,counter,m);
       
       arraysend[m]= array11[m+128*(counter)]+3*CMPLX(frame,0); 
       //printf("here is coredump self11");
       //printf("self11  real and imag is %f,%f\n",creal(arraysend[m]),cimag(arraysend[m]));
     }
    }

   // For self22 data
     if(floor(((counter/128)%3)==1))
    {
     for (int m=0;m<128;m=m+1)
     {
        //printf("\nm+128*counter %d , counter is %d , m is %d\n",m+128*counter,counter,m);
       
       arraysend[m]= array22[m+128*(counter-128)]+3*CMPLX(frame,0);
      //printf("here is coredump self22");

       //printf("self22  real and imag is %f,%f\n",creal(arraysend[m]),cimag(arraysend[m]));
     }
    }

      // For cross12 data
     if(floor(((counter/128)%3)==2))
    {
     for (int m=0;m<128;m=m+1)
     {
        //printf("\nm+128*counter %d , counter is %d , m is %d\n",m+128*counter,counter,m);
       
       arraysend[m]= array12[m+128*(counter-256)]+3*CMPLX(frame,0); 
        //printf("here is coredump cross12");
       //printf("cross real and imag is %f,%f\n",creal(arraysend[m]),cimag(arraysend[m]));
     }
    }
    


    



    //sleep(10);
    printf("\n======================================================================================================================================================================================\n");
    printf("frame number is %d\n",frame);
    for(int x=0;x<128;x++)
    {
	  
      //printf("sending array element %d is %f+j%f\n",x,creal(arraysend[x]),cimag(arraysend[x]));
      printf(" (%d) %f+j%f, ",x,creal(arraysend[x]),cimag(arraysend[x]));
      
    }     

    printf("\n=======================================================================================================================================================================================\n");
     z=sendto(sockfd, &arraysend, sizeof(arraysend), 0, (const struct sockaddr *) &servaddr,  sizeof(servaddr));
     //printf("\nZ value is %d\n",z);
     usleep(21000);
     if (z< 0) 
   {
    fprintf(stderr, "Could not send data\n");
   }
   else if (z=0)
   {
     printf("Nothing to send");
   }
   else if (z>0)
  {
    printf("Number of bytes sent is %d",z);
   }
   
   counter=counter+1;
 }  

  
 

  
    close(sockfd); 
    return 0; 
} 

