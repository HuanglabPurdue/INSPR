#include "cuda_runtime.h"
#include "definitions.h"
#include "kernel.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


__global__ void kernel_Localization(float *ParamIn, float *ParamNext, float *Convergence, float *FirstDev, float *SecondDev,
	int Nfit, int N_int, int FitBoxsize, float lambda, float SampleSpacingXY)
{
	const int tx = threadIdx.x;
	const int bx = blockIdx.x;
	const int BlockSize = blockDim.x;

	//Prevent read/write past end of array
	int j = BlockSize*bx + tx;
	if ((bx*BlockSize + tx) >= Nfit) return;

	
	float stepLimit[NP] = {0.03f, 0.03f, 0.06f, 400, 2}; // x,y,z step limits are in micron
	float x0_next[NPL];
	float dL_pos = 0, dL2_pos = 0;
	float dL_I, dL2_I; // photon and background
	float step[NPL];
	float rate = 1/(1 + lambda);
	float tmp;
	int s, p, k;
	// x,y,z
	for (p = 0; p < 3; p++)
	{
		for (s = 0; s < NCH; s++)		//Biplane; two channel; edited by FX
		{
			dL_pos += FirstDev[s*NP*Nfit + j*NP + p];
			dL2_pos += SecondDev[s*NP*Nfit + j*NP + p];
		}
		tmp = -1 * dL_pos / dL2_pos * rate;
		step[p] = fminf(fmaxf(tmp, -stepLimit[p]), stepLimit[p]);
	}

	for (s = 0; s < NCH; s++)
	{
		// photon
		dL_I = FirstDev[s*NP*Nfit + j*NP + 3];
		dL2_I = SecondDev[s*NP*Nfit + j*NP + 3];
		tmp = -1 * dL_I / dL2_I * rate;
		step[3 + s] = fminf(fmaxf(tmp, -stepLimit[3]), stepLimit[3]);
		// background
		dL_I = FirstDev[s*NP*Nfit + j*NP + 4];
		dL2_I = SecondDev[s*NP*Nfit + j*NP + 4];
		tmp = -1 * dL_I / dL2_I * rate;
		step[3 + NCH + s] = fminf(fmaxf(tmp, -stepLimit[4]), stepLimit[4]);		//??? step need change, edited by FX
	}

	x0_next[0] = ParamIn[NPL*j + 0] + step[0] * (-1 / SampleSpacingXY / N_int);
	x0_next[1] = ParamIn[NPL*j + 1] + step[1] * (-1 / SampleSpacingXY / N_int);
	for (k = 2; k < NPL; k++)
	{		
		x0_next[k] = ParamIn[NPL*j + k] + step[k];
	}

	for (s = 0; s < NCH; s++)
	{
		x0_next[3 + s] = (x0_next[3 + s] <= 100 ? 100 : x0_next[3 + s]); // intensity is not less than 100	
		x0_next[3 + NCH + s] = (x0_next[3 + NCH + s] <= 0 ? 0.01f : x0_next[3 + NCH + s]);// bg is not less than 0, edited by FX
	}
	x0_next[0] = fminf(fmaxf(x0_next[0], 4), FitBoxsize - 4);// xy shift is within fitting box
	x0_next[1] = fminf(fmaxf(x0_next[1], 4), FitBoxsize - 4);
	x0_next[2] = fminf(fmaxf(x0_next[2], -1.2), 1.2);//z position is within -1.4 to 1.4 um
	
	for (k = 0; k < NPL; k++) {
		ParamNext[NPL*j + k] = x0_next[k];
		Convergence[NPL*j + k] = x0_next[k] - ParamIn[NPL*j + k];
	}

}

__global__ void kernel_getdev(float *data, float *PSF, float *dPSFx, float *dPSFy, float *dPSFz, float *I, float *bg, int Nfit, int PSFsize,
	float *FirstDev, float *SecondDev)
{
	const int tx = threadIdx.x;
	const int bx = blockIdx.x;
	const int BlockSize = blockDim.x;
	
	//Prevent read/write past end of array
	int j = BlockSize*bx + tx;
	if ((bx*BlockSize + tx) >= Nfit) return;
	
	float dL[NP], dL2[NP];
	float psfI;
	int k, i;
	for (k = 0; k < NP; k++)
	{
		dL[k] = 0;
		dL2[k] = 0;
	}
	fundev(&data[j*PSFsize], &PSF[j*PSFsize], &dPSFx[j*PSFsize], I[j], I[j], bg[j], &dL[0], &dL2[0], PSFsize);
	fundev(&data[j*PSFsize], &PSF[j*PSFsize], &dPSFy[j*PSFsize], I[j], I[j], bg[j], &dL[1], &dL2[1], PSFsize);
	fundev(&data[j*PSFsize], &PSF[j*PSFsize], &dPSFz[j*PSFsize], I[j], I[j], bg[j], &dL[2], &dL2[2], PSFsize);
	fundev(&data[j*PSFsize], &PSF[j*PSFsize], &PSF[j*PSFsize], I[j], 1.0, bg[j], &dL[3], &dL2[3], PSFsize);
	for (int i = 0; i < PSFsize; i++)
	{
		psfI = PSF[j*PSFsize+i] * I[j] + bg[j];
		dL[4] += (data[j*PSFsize+i] / psfI - 1);
		dL2[4] += -1 * data[j*PSFsize+i] / psfI / psfI;

	}

	for (int k = 0; k < NP; k++)
	{
		FirstDev[NP * j + k] = dL[k];
		SecondDev[NP * j + k] = dL2[k];
	}
}

__device__ void fundev(float *data, float *psf, float *dpsf, float I, float Id, float bg, float *dL, float *dL2, int PSFsize)
{
	float psfI;
	for (int i = 0; i < PSFsize; i++)
	{
		psfI = psf[i] * I + bg;
		dL[0] += (data[i] / psfI - 1) * dpsf[i] * Id;
		dL2[0] += -1 * Id * Id * dpsf[i] * dpsf[i] * data[i] / psfI / psfI;
	}
}