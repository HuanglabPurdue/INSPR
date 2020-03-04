#include "cuda_runtime.h"
#include "definitions.h"
#include "kernel.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

__global__ void kernel_Err(float *PSF, float *Data, float *gainR, float *Error, int Nfit, int PSFSize){

	const int tx = threadIdx.x;
	const int bx = blockIdx.x;
	const int BlockSize = blockDim.x;
	//Prevent read/write past end of array
	int j = BlockSize*bx + tx;
	if ((bx*BlockSize + tx) >= Nfit) return;

	float sse = 0;
	float LLR = 0;
	for (int i = 0; i < PSFSize; i++)
	{
		sse += pow(PSF[j*PSFSize + i] - Data[j*PSFSize + i], 2);

		LLR += 2 * (PSF[j*PSFSize + i] - Data[j*PSFSize + i] -
			(Data[j*PSFSize + i] + gainR[j*PSFSize + i]) * log(PSF[j*PSFSize + i] + gainR[j*PSFSize + i])
			+ (Data[j*PSFSize + i] + gainR[j*PSFSize + i]) * log(Data[j*PSFSize + i] + gainR[j*PSFSize + i]));

		//mexPrintf("%.05f \n",toppsf[k*psfsize+i]);
	}
	Error[j * 2] += sse;
	Error[j * 2 + 1] += LLR;
}
