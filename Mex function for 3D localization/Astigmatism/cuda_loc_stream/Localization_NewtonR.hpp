#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cuda_runtime.h>
#include <math.h>
#include "definitions.h"
#include <mex.h>

//******************************************************
// extern C declarations of the kernels from wrapper.cu.  Yes you can put this in a
// header if you want.
extern "C" void kernel_Localization_wrapper(dim3, dim3, float *PSFa, float *PSFb, float *Data1, float *Data2,
	float *ParamIn, float *ParamNext, float * Convergence, int Nparam, int Nfit, int PSFSize, int FitBoxCenter,
	float deltax, float deltaz, float Iratio);

extern "C" void kernel_calCRLB_wrapper(dim3, dim3, float *PSFa, float *PSFb, float *ParamVar,
	int Nparam, int Nfit, int PSFSize, float delta, float Iratio);

class Localization_NewtonR{
public:

	float * PSFa;
	float * PSFb;
	float * Data1;
	float * Data2;
	float * ParamIn;
	float * ParamNext;
	float * ParamVar;
	float * Convergence;

	float * d_PSFa;
	float * d_PSFb;
	float * d_Data1;
	float * d_Data2;
	float * d_ParamIn;
	float * d_ParamNext;
	float * d_ParamVar;
	float * d_Convergence;

	int Nparam;
	int Nfit;
	int PSFSize;
	int Npsf;
	int FitBoxCenter;
	float Iratio;
	float Deltax;
	float Deltaz;

	void setup();
	void sendParams();
	void localization();
	void calCRLB();
	void cleanUp();
	void CudaError(const char *);


};