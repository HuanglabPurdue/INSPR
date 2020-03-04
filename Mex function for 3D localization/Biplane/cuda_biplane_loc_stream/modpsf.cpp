#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mex.h>
#include <cuda_runtime.h>
#include "definitions.h"
#include "kernel.h"

extern void kernel_PSF_image(dim3 numBlocks, dim3 threadsPerBlock,
							int Boxsize, float PixelSize,
							float *Xpos, float *Ypos, float *Z,
							float NA, float Lambda, 
							int *pCZ_n, int *pCZ_m,
							float *pCZ_real, float *pCZ_imag,
							int *pZernNum,
							int NzernS,
							float *rnm,
							float *pPSF,
							float n_med, float n_imm, float depth, float focalDis);

extern void kernel_convolve(dim3 numBlocks, dim3 threadsPerBlock, 
							int boxsize, int FOTFsize, 
							const float *OTFp,
							const float *oldPSF, 
							float * newPSF, 
							float *I, float *bg, float pixelsize);


void Modified_PSF(float *output, const float *x, const float *y, const float *z,
				int *d_zcn, int *d_zcm,
				float *d_zcr, float *d_zci,
				int *dZernNum,
				int NzernS,
				int numpsf, const float lambda, const float NA, const float pixelsize, const int boxsize,
				float *d_OTFparam,
				float *I, float *bg1, 
				float *d_rnm,
				float n_med, float n_imm, float depth, float focalDis)
{
	int i;
	int FOTFsize = 5;
	int OTFpsize = 4;
	int tt = (FOTFsize-1)/2;
	int boxsizeL = boxsize + FOTFsize - 1;
	const int sizeL = numpsf*boxsizeL*boxsizeL;
	const int size = numpsf*boxsize*boxsize;

	float *d_rawpsf=0;
	float *d_x=0;
	float *d_y=0;
	float *d_z=0;

	float *xtt=0;
	float *ytt=0; /* we need to add tt to both x and y to get the PSF centered properly in the larger grid */
	float *tempPSF=0; /* we also need a place to put the raw PSF image */

	xtt = (float *)mxMalloc(numpsf*sizeof(float));
	ytt = (float *)mxMalloc(numpsf*sizeof(float));
	tempPSF = (float *)mxMalloc(sizeL*sizeof(float));
	for(i=0;i<numpsf;i++)
	{
		xtt[i] = x[i] + tt;
		ytt[i] = y[i] + tt;
	}
	

	cudasafe(cudaMalloc((void**) &d_rawpsf, sizeL*sizeof(float)), "cudaMalloc d_PSF");
	cudasafe(cudaMemset(d_rawpsf, 0, sizeL*sizeof(float)), "set d_PSF to 0");

	cudasafe(cudaMalloc((void**) &d_x, numpsf*sizeof(float)), "cudaMalloc d_x");
	cudasafe(cudaMalloc((void**) &d_y, numpsf*sizeof(float)), "cudaMalloc d_y");
	cudasafe(cudaMalloc((void**) &d_z, numpsf*sizeof(float)), "cudaMalloc d_z");
	cudasafe(cudaMemcpy(d_x, xtt, numpsf*sizeof(float), cudaMemcpyHostToDevice), "copy x to device");
	cudasafe(cudaMemcpy(d_y, ytt, numpsf*sizeof(float), cudaMemcpyHostToDevice), "copy y to device");
	cudasafe(cudaMemcpy(d_z, z, numpsf*sizeof(float), cudaMemcpyHostToDevice), "copy z to device");

	dim3 numBlocks(numpsf, 1, 1);
	//mexPrintf("using %i blocks for raw psf calculation\n", numBlocks.x);
	dim3 threadsPerBlock(boxsizeL*boxsizeL, 1, 1);
	//mexPrintf("Block: %d, %d, %d\n", threadsPerBlock.x, threadsPerBlock.y, threadsPerBlock.z);

	kernel_PSF_image(numBlocks, threadsPerBlock,
					boxsizeL, pixelsize,
					d_y, d_x, d_z,
					NA, lambda, 
					d_zcn, d_zcm,
					d_zcr, d_zci,
					dZernNum,
					NzernS,
					d_rnm,
					d_rawpsf,
					n_med, n_imm, depth, focalDis);


	cudaError();
    cudaDeviceSynchronize(); //just for paranoia for now

	//cudasafe(cudaMemcpy(tempPSF, d_PSF, sizeL*sizeof(float), cudaMemcpyDeviceToHost), "copy PSF to host");
	//cudasafe(cudaFree(d_PSF), "free d_PSF");
	cudasafe(cudaFree(d_x), "free d_x");
	cudasafe(cudaFree(d_y), "free d_y");
	cudasafe(cudaFree(d_z), "free d_z");
	/* don't free d_rawpsf because we'll need it for the OTF calculation below */



	/*---------------------------------- end of raw PSF calculation--------------------------------------- */

	//float *d_rawpsf=0;
	float *d_PSF=0;
	float *d_I=0;
	float *d_bg1=0;

	cudasafe(cudaMalloc((void**) &d_PSF, size*sizeof(float)), "cudaMalloc d_PSF for OTF calculation");
	//cudasafe(cudaMalloc((void**) &d_rawpsf, sizeL*sizeof(float)), "cudaMalloc d_rawpsf");
	cudasafe(cudaMalloc((void**) &d_I, numpsf*sizeof(float)), "cudaMalloc d_I");
	cudasafe(cudaMalloc((void**) &d_bg1, numpsf*sizeof(float)), "cudaMalloc d_bg1");

	//cudasafe(cudaMemcpy(d_rawpsf, tempPSF, sizeL*sizeof(float), cudaMemcpyHostToDevice), "cudaMemcpy tempPSF to device");
	cudasafe(cudaMemcpy(d_I, I, numpsf*sizeof(float), cudaMemcpyHostToDevice), "cudaMemcpy I to device");
	cudasafe(cudaMemcpy(d_bg1, bg1, numpsf*sizeof(float), cudaMemcpyHostToDevice), "cudaMemcpy bg1 to device");

	dim3 numBlocks2(numpsf, 1, 1);
	//mexPrintf("using %i blocks for OTF adjustment\n", numBlocks2.x);
	dim3 threadsPerBlock2(boxsize*boxsize, 1, 1);
	//mexPrintf("OTF Block: %d, %d, %d\n", threadsPerBlock2.x, threadsPerBlock2.y, threadsPerBlock2.z);

	kernel_convolve(numBlocks2, threadsPerBlock2, 
					boxsize, FOTFsize, 
					d_OTFparam, 
					d_rawpsf, 
					d_PSF, 
					d_I, d_bg1, pixelsize);
	cudaError();
	cudaDeviceSynchronize(); //just for paranoia for now

	cudasafe(cudaMemcpy(output, d_PSF, size*sizeof(float), cudaMemcpyDeviceToHost), "copy final PSF to host");
	cudasafe(cudaFree(d_PSF), "free d_PSF");
	cudasafe(cudaFree(d_rawpsf), "free d_rawpsf");
	cudasafe(cudaFree(d_I), "free d_I");
	cudasafe(cudaFree(d_bg1), "free d_bg1");
	d_PSF = 0;
	d_rawpsf = 0;
	d_I = 0;
	d_bg1 = 0;

	//cudaDeviceReset();

	mxFree(xtt);
	mxFree(ytt);
	mxFree(tempPSF);
	xtt = ytt = tempPSF = 0;

}
