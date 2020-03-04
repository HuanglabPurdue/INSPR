#include "PSF_Interp2.hpp"

void PSF_Stream_Master::setup(){
	
	PSF = new PSF_Interp[NStreams];

	int NPSFperStream = floorf(NPSFs / NStreams);

	//Only want this once for all streams.  
	cudaMalloc(&d_SampledPSF, SampledPSF_XYPixels*SampledPSF_XYPixels*SampledPSF_ZPixels*sizeof(float));
	cudaMemcpyAsync(d_SampledPSF, SampledPSF, SampledPSF_XYPixels*SampledPSF_XYPixels*SampledPSF_ZPixels*sizeof(float), cudaMemcpyHostToDevice);

	//Loop through streams for setup
	for (int nn = 0; nn < NStreams; nn++){
		
		//Number of PSFs in this stream
		int NPSFthisStream = fminf(NPSFs - nn*NPSFperStream, NPSFperStream);

		//Propeties Stream independent
		PSF[nn].d_SampledPSF = d_SampledPSF;
		PSF[nn].SampledPSF_XYPixels = SampledPSF_XYPixels;
		PSF[nn].SampledPSF_ZPixels = SampledPSF_ZPixels;
		PSF[nn].SampleSpacingXY = SampleSpacingXY;
		PSF[nn].SampleSpacingZ = SampleSpacingZ;
		PSF[nn].StartX = StartX;
		PSF[nn].StartY = StartY;
		PSF[nn].StartZ = StartZ;
		PSF[nn].PSFs = PSFs;
		//Propeties Stream dependent
		PSF[nn].NPSFs = NPSFthisStream;

		//Allocate GPU memory and copy
		PSF[nn].setup();
	}

}

PSF_Stream_Master::~PSF_Stream_Master(){

	delete[] PSF;
}


void PSF_Stream_Master::calcPSF(float *X, float *Y, float *Z, float *I, float *Bg){

	int NPSFperStream = floorf(NPSFs / NStreams);
	
	for (int nn = 0; nn < NStreams; nn++){
		PSF[nn].NPSFs = NPSFperStream;
		PSF[nn].PSFs = &PSFs[nn*NPSFperStream * PSFSizeOut * PSFSizeOut];
		PSF[nn].sendPositions(&X[nn*NPSFperStream], &Y[nn*NPSFperStream], &Z[nn*NPSFperStream], &I[nn*NPSFperStream], &Bg[nn*NPSFperStream]);
	}
	for (int nn = 0; nn < NStreams; nn++){
		PSF[nn].PSFs = &PSFs[nn*NPSFperStream * PSFSizeOut * PSFSizeOut];
		PSF[nn].NPSFs = NPSFperStream;
		PSF[nn].calcPSF();
	}
	for (int nn = 0; nn < NStreams; nn++){
		PSF[nn].PSFs = &PSFs[nn*NPSFperStream * PSFSizeOut * PSFSizeOut];
		PSF[nn].NPSFs = NPSFperStream;
		PSF[nn].getPSF();
	}
	
}

void PSF_Stream_Master::cleanup(){
	for (int nn = 0; nn < NStreams; nn++)PSF[nn].cleanup();
	
}


PSF_Interp::~PSF_Interp(){

	//Any other clean-up

};

void PSF_Interp::setup(){

	//setup streams
	cudaStreamCreate(&CudaStream);

	//allocaate memory on device
	cudaMalloc(&d_PSFs, PSFSizeOut*PSFSizeOut*NPSFs*sizeof(float));
	cudaMalloc(&d_X, NPSFs*sizeof(float));
	cudaMalloc(&d_Y, NPSFs*sizeof(float));
	cudaMalloc(&d_Z, NPSFs*sizeof(float));
	cudaMalloc(&d_I, NPSFs*sizeof(float));
	cudaMalloc(&d_Bg, NPSFs*sizeof(float));
	CudaError("setup::Malloc");

};

void PSF_Interp::sendPositions(float * X, float * Y, float *Z, float *I, float *Bg){
	
	//copy from host to device
	cudaMemcpyAsync(d_X, X, NPSFs*sizeof(float), cudaMemcpyHostToDevice, CudaStream);
	cudaMemcpyAsync(d_Y, Y, NPSFs*sizeof(float), cudaMemcpyHostToDevice, CudaStream);
	cudaMemcpyAsync(d_Z, Z, NPSFs*sizeof(float), cudaMemcpyHostToDevice, CudaStream);
	cudaMemcpyAsync(d_I, I, NPSFs*sizeof(float), cudaMemcpyHostToDevice, CudaStream);
	cudaMemcpyAsync(d_Bg, Bg, NPSFs*sizeof(float), cudaMemcpyHostToDevice, CudaStream);
	CudaError("sendPositions::Memcpy");
}


void PSF_Interp::calcPSF(){

	//Each output pixel is the sum over 4x4 interpolated pixels.  
	//Interpolation is bilinear.  
	//Each block calculates one row of output pixels.   
	//This allows coallesed global memory reads accross threads.  
	int N_int = 4;		//There are N_int*N_int number of interpoloted points. 

	dim3 dimBlock(N_int,PSFSizeOut);
	dim3 dimGrid(PSFSizeOut, NPSFs);

	CudaError("calcPSF::setup Dims");
	mexPrintf("%d %d %d \n", N_int, PSFSizeOut, NPSFs);

	//kernel calls
	kernel_calcPSFPixel_wrapper(dimGrid, dimBlock, d_SampledPSF, SampledPSF_XYPixels, SampledPSF_XYPixels, d_PSFs, d_X, d_Y, d_Z, d_I, d_Bg,
		SampleSpacingXY, SampleSpacingZ, StartX, StartY, StartZ, N_int, PSFSizeOut, NPSFs, CudaStream);

	CudaError("calcPSF::kernel");

	

};

void PSF_Interp::getPSF(){

	cudaMemcpyAsync(PSFs, d_PSFs, NPSFs*PSFSizeOut*PSFSizeOut*sizeof(float), cudaMemcpyDeviceToHost, CudaStream);
	
	CudaError("calcPSF::memcpy");


}


void PSF_Interp::cleanup(){
	cudaFree(d_PSFs);
	cudaFree(d_X);
	cudaFree(d_Y);
	cudaFree(d_Z);
	cudaFree(d_I);
	cudaFree(d_Bg);
	//cudaStreamDestroy(CudaStream);

	CudaError("Cuda Free");

}


//*******************************************************************************************
void PSF_Interp::CudaError(const char * InStr) {
	/*!
	*  \brief A simple function to dump Cuda errors and abort the run.
	*/
	cudaError_t errornum;
	const char *str = 0;
	if (errornum = cudaGetLastError()) {
		str = cudaGetErrorString(errornum);
		cudaDeviceReset();
		mexErrMsgIdAndTxt("CudaTemplate:CUDA", "%s: %s\nYou should clear this function in MATLAB for proper operation.\n", InStr, str);
	}
}