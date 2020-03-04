#include "cuda_runtime.h"
#include "definitions.h"
#include "kernel.h"
#include "MatInvLib.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

__global__ void kernel_calCRLB(float *ParamF, float *ParamVar, int Nfit)
{
	const int tx = threadIdx.x;
	const int bx = blockIdx.x;
	const int BlockSize = blockDim.x;
	//Prevent read/write past end of array
	int j = BlockSize*bx + tx;
	if ((bx*BlockSize + tx) >= Nfit) return;

	int s, k;
	float FisherM[NPL*NPL];
	float LowerBi[NPL*NPL];
	float DiagLowerBi[NPL];

	for (k = 0; k < NPL*NPL; k++) FisherM[k] = 0;

	for (k = 0; k < NPL*NPL; k++)
	{
		for (s = 0; s < NCH; s++)	//Edited by FX
			FisherM[k] += ParamF[s*NPL*NPL*Nfit + j*NPL*NPL + k];
	}

	kernel_MatInvN(FisherM, LowerBi, DiagLowerBi, NPL);

	for (k = 0; k < NPL; k++) ParamVar[j*NPL + k] = DiagLowerBi[k];
}


__global__ void kernel_calFisherM(float *PSF, float *dPSFx, float *dPSFy, float *dPSFz, float *I, float *bg, float *gainR, float *ParamF,
	 int Q, int Nfit, int PSFSize){

	const int tx = threadIdx.x;
	const int bx = blockIdx.x;
	const int BlockSize = blockDim.x;
	//Prevent read/write past end of array
	int j = BlockSize*bx + tx;
	if ((bx*BlockSize + tx) >= Nfit) return;

	int t, k, i, s;
	float PSFa0;
	float funFi1[NPL];
	float tmp1;
	float FisherM[NPL*NPL];
	float w;

	for (i = 0; i < NPL*NPL; i++) FisherM[i] = 0;
	
		
	
	for (i = 0; i < PSFSize; i++)
	{   
		PSFa0 = PSF[j*PSFSize + i] * I[j] + bg[j] + gainR[j*PSFSize + i];
		//x 
		funFi1[0] = dPSFx[j*PSFSize + i] * I[j];
		//y
		funFi1[1] = dPSFy[j*PSFSize + i] * I[j];
		//z
		funFi1[2] = dPSFz[j*PSFSize + i] * I[j];

		for (s = 0; s < NCH; s++)	//edited by FX
		{
			w = (Q == s ? 1.0 : 0.0);
			// I 
			funFi1[3 + s] = PSF[j*PSFSize + i] * w;
			// bg
			funFi1[3 + NCH + s] = w;
		}
		
		for (t = 0; t < NPL; t++)
		{
			for (k = 0; k < NPL; k++)
			{
				tmp1 = funFi1[t] * funFi1[k] / fmaxf(PSFa0, 1e-4f);
				FisherM[t*NPL + k] += tmp1;
			}
		}
		
	}

	for (k = 0; k < NPL*NPL; k++) ParamF[j*NPL*NPL + k] = FisherM[k];
}

