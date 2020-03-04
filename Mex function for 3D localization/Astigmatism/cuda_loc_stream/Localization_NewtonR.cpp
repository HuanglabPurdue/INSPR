#include "Localization_NewtonR.hpp"


void Localization_NewtonR::setup(){

	//This resets the card and clears any remaining error codes. We should not do this if we eventually make multiple calls sending in sampled stack only
	// once since that will clear memory. 
	//cudaDeviceReset();

	//allocate memory on device
	cudaMalloc(&d_PSFa, Npsf*Nfit*PSFSize*sizeof(float));
	cudaMalloc(&d_PSFb, Npsf*Nfit*PSFSize*sizeof(float));
	cudaMalloc(&d_Data1, Nfit*PSFSize*sizeof(float));
	cudaMalloc(&d_Data2, Nfit*PSFSize*sizeof(float));
	cudaMalloc(&d_ParamIn, Nfit*Nparam*sizeof(float));
	cudaMalloc(&d_ParamVar, Nfit*Nparam*sizeof(float));
	cudaMalloc(&d_ParamNext, Nfit*Nparam*sizeof(float));
	cudaMalloc(&d_Convergence, Nfit*Nparam*sizeof(float));

	CudaError("Malloc");

	//copy from host to device
	cudaMemcpy(d_Data1, Data1, Nfit*PSFSize*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_Data2, Data2, Nfit*PSFSize*sizeof(float), cudaMemcpyHostToDevice);
	CudaError("Memcpy");

};

void Localization_NewtonR::sendParams(){
	//copy from host to device
	cudaMemcpy(d_PSFa, PSFa, Npsf*Nfit*PSFSize*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_PSFb, PSFb, Npsf*Nfit*PSFSize*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_ParamIn, ParamIn, Nfit*Nparam*sizeof(float), cudaMemcpyHostToDevice);
	CudaError("Memcpy");
}
void Localization_NewtonR::cleanUp(){

	//free all memory we allocated. 
	cudaFree(d_PSFa);
	cudaFree(d_PSFb);
	cudaFree(d_Data1);
	cudaFree(d_Data2);
	cudaFree(d_ParamIn);
	cudaFree(d_ParamNext);
	cudaFree(d_ParamVar);
	cudaFree(d_Convergence);
	CudaError("Cuda Free");

};

void Localization_NewtonR::localization(){
	//Max threads per block is 1024 = 32^2.  Since we might want bigger than 32*32, lets have each block do one row of pixels. 
	int BlockSize = BSZ;
	int NBlock = Nfit / BlockSize + 1;
	dim3 dimBlock(BlockSize);
	dim3 dimGrid(NBlock);

	//kernel calls
	kernel_Localization_wrapper(dimGrid, dimBlock, d_PSFa, d_PSFb, d_Data1, d_Data2,
		d_ParamIn, d_ParamNext, d_Convergence, Nparam, Nfit, PSFSize, FitBoxCenter,
		Deltax, Deltaz, Iratio);

	cudaMemcpy(ParamNext, d_ParamNext, Nfit*Nparam*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(Convergence, d_Convergence, Nfit*Nparam*sizeof(float), cudaMemcpyDeviceToHost);

	CudaError("kernel_Localization_wrapper");

}

void Localization_NewtonR::calCRLB(){
	//Max threads per block is 1024 = 32^2.  Since we might want bigger than 32*32, lets have each block do one row of pixels. 
	int BlockSize = BSZ;
	int NBlock = Nfit / BlockSize + 1;
	dim3 dimBlock(BlockSize);
	dim3 dimGrid(NBlock);

	//kernel calls
	kernel_calCRLB_wrapper(dimGrid, dimBlock, d_PSFa, d_PSFb, d_ParamVar,
		Nparam, Nfit, PSFSize, Deltax, Iratio);

	cudaMemcpy(ParamVar, d_ParamVar, Nfit*Nparam*sizeof(float), cudaMemcpyDeviceToHost);
	CudaError("kernel_Localization_wrapper");

}

void Localization_NewtonR::CudaError(const char * InStr) {
	/*!
	*  \brief A simple function to dump Cuda errors and abort the run.
	*/
	cudaError_t errornum;
	const char *str = 0;
	if (errornum = cudaGetLastError()) {
		str = cudaGetErrorString(errornum);
		mexErrMsgIdAndTxt("CudaTemplate:CUDA", "%s: %s\nYou should clear this function in MATLAB for proper operation.\n", InStr, str);
	}
}