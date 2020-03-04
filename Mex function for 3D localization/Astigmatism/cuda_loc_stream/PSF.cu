#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mex.h>
#include <cuda_runtime.h>
#include "kernel.h"
#include "definitions.h"


__device__ void fkr(float factor, int n, int m, 
				int kk, float k, float krmax, 
				float rho, float phi, 
				float Z,
				float * rnm, 
				float pCZ_real, float pCZ_imag, 
				int ZernNum, 
				float *IntRe, float *IntIm,
				float Lambda,float n_med,float n_imm,
				float depth,float focalDis,float kccd)
{	
	float ipower[] = {0, 0, 0, 0};
	float	dkr=krmax/SPD;
	float   kr=kk*dkr;
	//float	argument = 2*PI*kr*rho; // argument of the besselfunctions inside the Zernike sum
	float focusfactor, focusCCD,aberPhase,aberPhaseRe,aberPhaseIm;  // argument of the "fourier transform" part of the integral
	float Cos_theta1,Cos_theta3;
	float	fftRe, fftIm;
	float	sumRe,sumIm;
	
	int		rnmbase=(n*n+n)/2+(n-m);
	
	int		rnmcount=rnmbase*SPD+kk;
	int Qcos=(ZernNum==(n*n+2*(n-m)) ? 1 : 0);
	int Q=1-pow(kr*Lambda/n_med,2)>0 ? 1:0;
	
	switch (Q)
		{
			case 1:
				Cos_theta1=sqrtf(1-pow(kr*Lambda/n_med,2));
				Cos_theta3=sqrtf(1-pow(kr*Lambda/n_imm,2));
				aberPhase=2*PI/Lambda*(depth*n_med*Cos_theta1-depth*n_imm*n_imm/n_med*Cos_theta3);
				focusfactor=2*PI*k*Z*Cos_theta1;
				focusCCD=2*PI*kccd*focalDis*Cos_theta3;
				fftRe = dkr*kr*cosf(focusfactor+aberPhase+focusCCD);
				fftIm = dkr*kr*sinf(focusfactor+aberPhase+focusCCD); //Do not Fudge this negative sign!!!
				break;
			case 0:
				Cos_theta1=sqrtf(pow(kr*Lambda/n_med,2)-1); // imaginary 
				Cos_theta3=sqrtf(1-pow(kr*Lambda/n_imm,2));
				aberPhaseIm=2*PI*depth*k*Cos_theta1;
				aberPhaseRe=-2*PI*depth*k*pow(n_imm/n_med,2)*Cos_theta3;
				focusfactor=2*PI*k*Z*Cos_theta1;
				focusCCD=2*PI*kccd*focalDis*Cos_theta3;
				fftRe = dkr*kr*exp(-aberPhaseIm-focusfactor)*cosf(aberPhaseRe+focusCCD);
				fftIm = dkr*kr*exp(-aberPhaseIm-focusfactor)*sinf(aberPhaseRe+focusCCD); //Do not Fudge this negative sign!!!
				break;
		}

		factor*=rnm[rnmcount];	
	
		switch (Qcos){
		case 1:

			ipower[m % 4] += factor*pCZ_real*cosf(m*phi);
			ipower[(m+1) % 4] += factor*pCZ_imag*cosf(m*phi);
			break;
		case 0:
			ipower[m % 4] += factor*pCZ_real*sinf(m*phi);
			ipower[(m+1) % 4] += factor*pCZ_imag*sinf(m*phi);
			break;
	}
			
	sumRe = ipower[0] - ipower[2];
	sumIm = ipower[1] - ipower[3];

	IntRe[0] = fftRe*sumRe - fftIm*sumIm;
	IntIm[0] = fftIm*sumRe + fftRe*sumIm;	

		// the above places the +1, +i, -1, and -i terms (in that order) in separate columns to be added (subtracted) later
		// note that I am not accounting for the negatives resulting from powers of i--this will be taken care of once we exit the loop
}

__device__ float PSF_point(float rho, float phi,
						   float Z,
						   float NA,float Lambda,
						   int *pCZ_n, int *pCZ_m, 
						   float *pCZ_real,float *pCZ_imag, 
						   int *pZernNum, 
						   float *rnm,  
						   int NzernS,
						   float n_med,float n_imm,float depth,float focalDis)
{
	float integralRe = 0, integralIm = 0;	
	float prefactor = 2.0f*PI;
	int	 nn, mm, kk, ii ;
	float IntRe,IntIm;
	float CZ_real,CZ_imag;
	float factor;
	float argument;
	float dkr,kr;
	int ZernN;
	float k= n_med/Lambda; 
    float krmax= NA/Lambda; //krRadPoly is defined as kr/krmax and is used for the radial polynomial which must be normalized...
	float kccd=n_imm/Lambda;
	
	dkr=krmax/SPD;

	for(ii=0;ii<NzernS;ii++){
		mm=pCZ_m[ii];
		nn=pCZ_n[ii];
		CZ_real=pCZ_real[ii];
		CZ_imag=pCZ_imag[ii];
		ZernN=pZernNum[ii];
		switch (mm)
		{
		case 0:
			//precalculate kk==0 case
			kr = 0;
			// mm == 0
			factor = 1; //besselj(0,0)

					fkr(factor, nn, 0, 
						0, k, krmax, 
						rho, phi, 
						Z, 
						rnm, 
						CZ_real, CZ_imag,
						ZernN, 
						&IntRe, &IntIm,
						Lambda,n_med,n_imm,depth,focalDis,kccd);

			integralRe += IntRe;
			integralIm += IntIm;	
			for (kk=1;kk<SPD;kk++){
				kr=kk*dkr;
				argument = 2*PI*kr*rho;
				// calculate mm==0 case
				factor = j0f(argument);

					fkr(factor, nn, 0, 
						kk, k, krmax, 
						rho, phi, 
						Z, 
						rnm, 
						CZ_real, CZ_imag,
						ZernN, 
						&IntRe, &IntIm,
						Lambda,n_med,n_imm,depth,focalDis,kccd);

				integralRe += IntRe;
				integralIm += IntIm;	
			}
			break;
		case 1:
			// mm == 1
			kr=0;
			factor = 0; //besselj(1,0)

					fkr(factor, nn, 1, 
						0, k, krmax, 
						rho, phi, 
						Z, 
						rnm, 
						CZ_real, CZ_imag,
						ZernN, 
						&IntRe, &IntIm,
						Lambda,n_med,n_imm,depth,focalDis,kccd);

			integralRe += IntRe;
			integralIm += IntIm;	

			for (kk=1;kk<SPD;kk++){
				kr=kk*dkr;
				argument = 2*PI*kr*rho;
				factor = j1f(argument);

					fkr(factor, nn, 1, 
						kk, k, krmax, 
						rho, phi, 
						Z, 
						rnm, 
						CZ_real, CZ_imag,
						ZernN, 
						&IntRe, &IntIm,
						Lambda,n_med,n_imm,depth,focalDis,kccd);

				integralRe += IntRe;
				integralIm += IntIm;	
			}
			break;
		default:// mm > 1

			factor = 0; //besselj(n,0)

					fkr(factor, nn, mm, 
						0, k, krmax,
						rho, phi, 
						Z, 
						rnm, 
						CZ_real, CZ_imag,
						ZernN, 
						&IntRe, &IntIm,
						Lambda,n_med,n_imm,depth,focalDis,kccd);

			integralRe += IntRe;
			integralIm += IntIm;	
			for (kk=1;kk<SPD;kk++){
				kr=kk*dkr;
				argument = 2*PI*kr*rho;
				factor = jnf(mm,argument);

					fkr(factor, nn, mm, 
						kk, k, krmax, 
						rho, phi, 
						Z, 
						rnm, 
						CZ_real, CZ_imag,
						ZernN, 
						&IntRe, &IntIm,
						Lambda,n_med,n_imm,depth,focalDis,kccd);

				integralRe += IntRe;
				integralIm += IntIm;	

			}
			break;
		}

	}

	integralRe *= prefactor;
	integralIm *= prefactor;

	return integralRe*integralRe + integralIm*integralIm;
}

__device__ float OTFconvolve(const float * oldpsf, 
							 int boxsizeL, int tt, 
							 const float * OTFp, 
							 int x, int y, 
							 float I, float bg, float eps)
{
	int i,j;
	float osfpix=0.0; /* this is the return value */
	float expfactor1 = -2*OTFp[1]*OTFp[1]*PI*PI*eps*eps;
	float expfactor2 = -2*OTFp[2]*OTFp[2]*PI*PI*eps*eps;

	for (i=-tt;i<=tt;i++)
		for(j=-tt;j<=tt;j++)
			osfpix += 2*PI*OTFp[0]*OTFp[1]*OTFp[2]*eps*eps*expf(i*i*expfactor1)*expf(j*j*expfactor2)*oldpsf[((y+tt)-i)*boxsizeL+(x+tt)-j];

	return osfpix*I+bg;
}

/* Note that tempPSF is a global array that needs to be big enough to hold a PSF measuring
	 boxsizeL*boxsizeL*sizeof(float), whereas pPSF, which is the corrected PSF to be returned
	 is smaller, measuring boxsize*boxsize*sizeof(float) */
__global__ void kernel_PSF_image(int Boxsize,float PixelSize,
								 float *Xpos,float *Ypos,float *Z,
								 float NA,float Lambda,
								 int *pCZ_n, int *pCZ_m, 
								 float *pCZ_real,float *pCZ_imag,
								 int *pZernNum, 
								 int NzernS, 
								 float *rnm, 
								 float *pPSF,
								 float n_med,float n_imm,float depth,float focalDis)
{
	float X,Y,Rho,Phi;
	float norm,FTsize;//,maxK,dK;
	//int SamplesPerDim;
	
	//SamplesPerDim=SPD;
	//maxK=NA/Lambda;
	//dK=2*maxK/SamplesPerDim;
	
	// normalization
	//FTsize=1/PixelSize/dK; //FT equivalent size
	//norm=FTsize*FTsize*PI*(SamplesPerDim/2)*(SamplesPerDim/2);
	FTsize=PixelSize;
	norm=FTsize*FTsize*FTsize*FTsize;

	const int pixel = blockDim.x*threadIdx.y + threadIdx.x;
	const int pixelx = pixel % Boxsize;
	const int pixely = pixel / Boxsize;
	//const int pixel = blockIdx.x * blockDim.x*blockDim.x + blockIdx.y*blockDim.x+threadIdx.x;
	//const int pixelx = threadIdx.x;
	//const int pixely = blockIdx.y;

	X=-((float)pixelx-Xpos[blockIdx.x])*PixelSize;
	Y=-((float)pixely-Ypos[blockIdx.x])*PixelSize;
	Rho = sqrtf(X*X + Y*Y);
	Phi = atan2f(X, Y);
	pPSF[blockIdx.x*blockDim.x + pixel]=PSF_point(Rho,Phi,
						  Z[blockIdx.x],
						  NA,Lambda,
						  pCZ_n,pCZ_m,
						  pCZ_real,pCZ_imag,
						  pZernNum,
						  rnm,
						  NzernS,
						  n_med,n_imm,depth,focalDis)*norm;
}

__global__ void kernel_convolve(int boxsize, int FOTFsize, 
								const float * OTFp, 
								const float * oldPSF, 
								float * newPSF, 
								float *I, float *bg, float pixelsize)
{
				int boxsizeL = boxsize + FOTFsize - 1;
				int tt = (FOTFsize-1)/2;

				const int pixel = blockDim.x*threadIdx.y + threadIdx.x; //looks like you're only doing 1d here
				const int pixelx = pixel % boxsize;
				const int pixely = pixel / boxsize;
				//const int pixel = blockIdx.x * blockDim.x*blockDim.x + blockIdx.y*blockDim.x+threadIdx.x;

				//const int pixelx = threadIdx.x;
				//const int pixely = blockIdx.y;

				newPSF[blockIdx.x*blockDim.x + pixel] = OTFconvolve(&oldPSF[blockIdx.x*boxsizeL*boxsizeL], 
											boxsizeL, tt, 
											OTFp, 
											pixelx, pixely, 
											I[blockIdx.x], bg[blockIdx.x], pixelsize);
}
