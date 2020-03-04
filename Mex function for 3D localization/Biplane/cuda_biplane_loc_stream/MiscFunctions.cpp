#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mex.h>
#include <cuda_runtime.h>
#include "definitions.h"


void LU_decom(float *LU, int *mutate, float *A, int sz);

void GetInitialGuess(int Start, int boxsize, int channelNum,int s, float Iratio,  
	float *channel1, float *channel2, float *bg1, float *bg2, float *p, float *x0)
{
	int i, j,k;
	float comx1=0,comy1=0,area=0;
	float BG,BG1,BG2;//BG;
	float Ip1=0,Ip2=0,idefo,intensity1=0;
	int psfsize=boxsize*boxsize;
	int x0Size=6;
	// calculate center of mass
	for(k=0;k<channelNum;k++)
	{
		area=0;
		comx1=0;
		comy1=0;
		for(i=0;i<4;i++)
		{
			bg1[i]=0;
			bg2[i]=0;
		}
		Ip1=0;
		Ip2=0;
		intensity1=0;
		for(i=Start-2;i<=Start+2;i++)
			for(j=Start-2;j<=Start+2;j++)
			{
				area+=channel1[k*psfsize+i*boxsize+j];
				comx1+=channel1[k*psfsize+i*boxsize+j]*i;
				comy1+=channel1[k*psfsize+i*boxsize+j]*j;
			}
			comx1=comx1/area;
			comy1=comy1/area;

			// calculate bg1 , bg2
			for(i=0;i<boxsize;i++)
			{
				bg1[0]+=channel1[k*psfsize+i];
				bg1[1]+=channel1[k*psfsize+(boxsize-1)*boxsize+i];
				bg1[2]+=channel1[k*psfsize+i*boxsize];
				bg1[3]+=channel1[k*psfsize+i*boxsize+boxsize-1];
				bg2[0]+=channel2[k*psfsize+i];
				bg2[1]+=channel2[k*psfsize+(boxsize-1)*boxsize+i];
				bg2[2]+=channel2[k*psfsize+i*boxsize];
				bg2[3]+=channel2[k*psfsize+i*boxsize+boxsize-1];
			}

			for(i=0;i<4;i++)
			{
				bg1[i]=abs(bg1[i]/boxsize);
				bg2[i]=abs(bg2[i]/boxsize);
			}

			BG1=bg1[0];
			BG2=bg2[0];
			for(i=1;i<4;i++)
			{
				BG1=min(BG1,bg1[i]);
				BG2=min(BG2,bg2[i]);
			}
			BG=min(BG1,BG2);
			//BG=min(BG1,BG2); // This calculation not necessary?

			// estimate z
			for(i=Start-s;i<=Start+s;i++)
				for(j=Start-s;j<=Start+s;j++)
				{
					Ip1+=channel1[k*psfsize+i*boxsize+j]-BG;
					Ip2+=channel2[k*psfsize+i*boxsize+j]-BG;
				}

				idefo=Ip2/Ip1*p[0]+p[1];
				idefo=min(max(idefo,-1),1);;

				// estimate intensity1
				for(i=0;i<psfsize;i++) 
					intensity1+=channel1[k*psfsize+i]-BG1+(channel2[k*psfsize+i]-BG2)/Iratio;
				intensity1=intensity1/2;

				x0[k*x0Size+0]=comx1;
				x0[k*x0Size+1]=comy1;
				x0[k*x0Size+2]=intensity1;
				x0[k*x0Size+3]=BG1;
				x0[k*x0Size+4]=idefo;
				x0[k*x0Size+5]=BG2;
	}
}

void MatrixInv(float *M_inv,float *A, int sz)
{
	int i,j,k;
	float tmp1;
	float *LU=0;
	int *mutate=0;
	float *Y=0;
	float *Identity=0;


	LU=(float*)malloc(sz*sz*sizeof(float));
	mutate=(int*)malloc(sz*sizeof(int));
	Y=(float*)malloc(sz*sz*sizeof(float));
	Identity=(float*)malloc(sz*sz*sizeof(float));

	for(i=0;i<sz;i++)
	{
		for(j=0;j<sz;j++)
		{
			if(i==j)
				Identity[i*sz+j]=1;
			else
				Identity[i*sz+j]=0;
		}
	}

	
    LU_decom(LU,mutate,A,sz);
	

	for(j=0;j<sz;j++)
	{
		for(i=0;i<sz;i++)
		{
			tmp1=0;
			if(i==0)
			{
				Y[i*sz+j]=Identity[mutate[i]*sz+j];
			}
			else
			{
				for(k=0;k<=i-1;k++) tmp1+=LU[i*sz+k]*Y[k*sz+j];
				Y[i*sz+j]=Identity[mutate[i]*sz+j]-tmp1;
			}
		}

		for(i=sz-1;i>=0;i--)
		{
			
			tmp1=0;
			if(i==(sz-1))
			{
				M_inv[i*sz+j]=Y[i*sz+j]/LU[i*sz+i];
			}
			else
			{
				for(k=i+1;k<sz;k++) tmp1+=LU[i*sz+k]*M_inv[k*sz+j];
				M_inv[i*sz+j]=(Y[i*sz+j]-tmp1)/LU[i*sz+i];
				//mexPrintf("LU : %0.5f \n",LU[i*sz+i]);
			}
		}
	}

	

	free(LU);
	free(mutate);
	free(Y);
	free(Identity);
	LU = 0;
	mutate = 0;
	Y=0;
	Identity=0;

}


void LU_decom(float *LU, int *mutate, float *A, int sz)
{
	int i,j,k,s;
	int n,tt;
	float tmp1,p;
	
	float *row0=0;

	row0=(float*)malloc(sz*sizeof(float));


	for(s=0;s<sz;s++) mutate[s]=s;

	for(j=0;j<sz;j++)
	{
		for(i=0;i<=j;i++)
		{
			tmp1=0;
			if(i==0)
			{
				LU[i*sz+j]=A[i*sz+j];
			}
			else
			{
				for(k=0;k<=i-1;k++) tmp1+=LU[i*sz+k]*LU[k*sz+j];
					LU[i*sz+j]=A[i*sz+j]-tmp1;
					//mexPrintf("%0.5f \n",LU[i*sz+j]);
			}
			
			if(i==j)
			{
				p=abs(LU[i*sz+j]);
				n=j;
			}
		}

		for(i=j+1;i<sz;i++)
		{
			tmp1=0;
			if(j==0)
			{
				LU[i*sz+j]=A[i*sz+j];
			}
			else
			{
				for(k=0;k<=j-1;k++) tmp1+=LU[i*sz+k]*LU[k*sz+j];
				LU[i*sz+j]=A[i*sz+j]-tmp1;
			}

			if(abs(LU[i*sz+j])>p)
			{
				p=abs(LU[i*sz+j]);
				n=i;
			}

			//if(p==0)
			if(n!=j)
			{
				for(s=0;s<sz;s++) row0[s]=LU[j*sz+s];
				for(s=0;s<sz;s++) LU[j*sz+s]=LU[n*sz+s];
				for(s=0;s<sz;s++) LU[n*sz+s]=row0[s];

				for(s=0;s<sz;s++) row0[s]=A[j*sz+s];
				for(s=0;s<sz;s++) A[j*sz+s]=A[n*sz+s];
				for(s=0;s<sz;s++) A[n*sz+s]=row0[s];

				tt=mutate[j];
				mutate[j]=mutate[n];
				mutate[n]=tt;
				n=j;
			}

		}
		//tmp2=LU[j*sz+j];
		//if (tmp2==0) LU[j*sz+j]=1e-6;
		
		LU[j*sz+j]=max(abs(LU[j*sz+j]),1e-5f)*LU[j*sz+j]/abs(LU[j*sz+j]);
		for(i=j+1;i<sz;i++) {LU[i*sz+j]/=LU[j*sz+j];
		                      
		}


		
	}

	//for(j=0;j<sz;j++)
	//mexPrintf("%0.5f \n",LU[j*sz+j]);

	free(row0);
	row0=0;
}
