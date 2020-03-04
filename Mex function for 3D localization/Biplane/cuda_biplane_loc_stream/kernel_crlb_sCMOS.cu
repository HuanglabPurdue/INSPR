#include "cuda_runtime.h"
#include "definitions.h"
#include "kernel.h"
#include "MatInvLib.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


__global__ void kernel_calCRLB_sCMOS(float *PSF, float *GainR, float *ParamVar,
	int Nparam, int Nfit, int PSFSize, float delta){

	const int tx = threadIdx.x;
	const int bx = blockIdx.x;
	const int BlockSize = blockDim.x;
	//Prevent read/write past end of array
	int j = BlockSize*bx + tx;
	if ((bx*BlockSize + tx) >= Nfit) return;

	int t, k, i;
	float PSFa0;
	float funFi1[NP];
	float tmp1;
	float FisherM[NP*NP];
	float LowerBi[NP*NP];
	float DiagLowerBi[NP];

	for (i = 0; i < Nparam*Nparam; i++) FisherM[i] = 0;


	for (i = 0; i < PSFSize; i++)
	{
		PSFa0 = PSF[Vadex(j, 0)];
		//x 
		funFi1[0] = (PSF[Vadex(j, 2)] - PSFa0) / delta;
		//y
		funFi1[1] = (PSF[Vadex(j, 3)] - PSFa0) / delta;
		//I plane1
		funFi1[2] = PSF[Vadex(j, 1)];
		// bg1
		funFi1[3] = 1;
		//z
		funFi1[4] = (PSF[Vadex(j, 4)] - PSFa0) / delta;

		for (t = 0; t < Nparam; t++)
		{
			for (k = 0; k < Nparam; k++)
			{
				tmp1 = funFi1[t] * funFi1[k] / fmaxf(PSFa0 + GainR[PSFSize*j + i], 1e-4f);// add variance-gain ratio: v/g^2
				FisherM[t*Nparam + k] += tmp1;
			}
		}
	}
	kernel_MatInvN(FisherM, LowerBi, DiagLowerBi, Nparam);
	for (k = 0; k < Nparam; k++) ParamVar[j*Nparam + k] = DiagLowerBi[k];
}