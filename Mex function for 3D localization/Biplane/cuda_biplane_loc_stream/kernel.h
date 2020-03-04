
#include "Model_Localization.hpp"

#ifndef KERNEL_H
#define KERNEL_H


__global__ void kernel_calcPSFPixel(float *SamplePSF, float *SamplePSFx, float *SamplePSFy, float *SamplePSFz, float *SamplePSFxy, float *SamplePSFxz, float *SamplePSFyz, float *SamplePSFxyz,
	int SizeX, int SizeY, float *PSFs, float *dPSFx, float *dPSFy, float *dPSFz, float * X, float *Y, float *Z,
	float SampleSpacingXY, float SampleSpacingZ, float StartX, float StartY, float StartZ, int N_int, int PSFSizeOut, int NPSFs, float *d_tMat, float *d_offset_int, int flag);	//spline parameters, affine parameters

__global__ void kernel_calcPSFI(float *psf, float *psfI, float *I, float *bg, int Nfit, int PSFsize);

__global__ void kernel_Localization(float *ParamIn, float *ParamNext, float *Convergence, float *FirstDev, float *SecondDev,
	int Nfit, int N_int, int FitBoxCenter, float lambda, float SampleSpacingXY);

__global__ void kernel_getdev(float *data, float *PSF, float *dPSFx, float *dPSFy, float *dPSFz, float *I, float *bg, int Nfit, int PSFsize,
	float *FirstDev, float *SecondDev);

__global__ void kernel_calCRLB(float *ParamF, float *ParamVar, int Nfit);
	
__global__ void kernel_calFisherM(float *PSF, float *dPSFx, float *dPSFy, float *dPSFz, float *I, float *bg, float *ParamF,
	int Q, int Nfit, int PSFSize, float *d_tMat);	//affine

__global__ void kernel_Err(float *PSF, float *Data, float *Error, int Nfit, int PSFSize);


__device__ void gencoeff(float f1, float f2, float df1, float df2, float dx, float x1, float x2, float *a, float *b, float *c, float *d);

__device__ float evalspline(float a, float b, float c, float d, float x);

__device__ void dev_spline(float a, float b, float c, float x, float *df);

__device__ void getsampleVal(float *SamplePSF, int ind, int SizeX, float *F);

__device__ void fundev(float *data, float *psf, float *dpsf, float I, float Id, float bg, float *dL, float *dL2, int PSFsize);

__device__ void interpXY(float *SamplePSF, float *SamplePSFx, float *SamplePSFy, float *SamplePSFxy, int tmp2, int SizeX, float SampleSpacingXY, float X1x, float X1y,
	float X_thread, float Y_thread, float *F2d, float *dF2d_y, float *dF2d_y2, float *dF2d_x, float *dF2d_x2);


#endif