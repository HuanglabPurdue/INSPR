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
extern "C" void kernel_calcPSFPixel_wrapper(dim3, dim3, float *SampledPSF, int SizeX, int SizeY, float *PSFs, float * X, float *Y, float *Z, float *I, float *Bg,
	float SampleSpacingXY, float SampleSpacingZ, float StartX, float StartY, float StartZ, int N_int, int PSFSizeOut, int NPSFs, cudaStream_t CudaStream);




class PSF_Interp{
public:

	//Propeties 
	float * PSFs;

	float * d_SampledPSF;
	float * d_PSFs;

	float * d_X;
	float * d_Y;
	float * d_Z;
	float * d_I;
	float * d_Bg;

	int NPSFs;

	int SampledPSF_XYPixels;  //Size of sampled PSF
	int SampledPSF_ZPixels;

	float SampleSpacingXY;
	float SampleSpacingZ;
	float StartX;
	float StartY;
	float StartZ;

	int PSFSizeOut=16;
	
	//Cuda Properties
	cudaStream_t  CudaStream;
	int CudaDevice=0; 

	//Methods
	void setup();
	void sendPositions(float *, float *, float *, float *, float *);
	void calcPSF();
	void getPSF();
	void cleanup();

	~PSF_Interp();
	void CudaError(const char *);

	
};

class PSF_Stream_Master{

public:
	int NStreams = 2;
	float * SampledPSF;
	float * PSFs;
	float * d_SampledPSF;
	float * d_PSFs;
	float * X;
	float * Y;
	float * Z;
	float * I;
	float * Bg;

	int SampledPSF_XYPixels;  //Size of sampled PSF
	int SampledPSF_ZPixels;

	float SampleSpacingXY;
	float SampleSpacingZ;
	float StartX;
	float StartY;
	float StartZ;
	int NPSFs;
	int PSFSizeOut = 16;

	PSF_Interp *PSF;

	void setup();
	void calcPSF(float *, float *, float *, float *, float *);
	void cleanup();
	~PSF_Stream_Master();

};
