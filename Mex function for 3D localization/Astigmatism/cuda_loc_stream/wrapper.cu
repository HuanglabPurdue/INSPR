

#include "definitions.h"
#include "kernel.h"



//*******************************************************************************************
extern "C" void kernel_calCRLB_wrapper(dim3 dimGrid, dim3 dimBlock, float *ParamF, float *ParamVar, int Nfit)	
{
	kernel_calCRLB << <dimGrid, dimBlock >> >(ParamF, ParamVar, Nfit);
		
}

extern "C" void kernel_calFisherM_wrapper(dim3 dimGrid, dim3 dimBlock, float *PSF, float *dPSFx, float *dPSFy, float *dPSFz, float *I, float *bg, float *gainR, float *ParamF,
	int Q, int Nfit, int PSFSize)
{
	kernel_calFisherM << <dimGrid, dimBlock >> >(PSF, dPSFx, dPSFy, dPSFz, I, bg, gainR, ParamF, Q, Nfit, PSFSize);
}

extern "C" void kernel_Localization_wrapper(dim3 dimGrid, dim3 dimBlock, float *ParamIn, float *ParamNext, float *Convergence, float *FirstDev, float *SecondDev,
	int Nfit, int N_int, int FitBoxCenter, float Lambda, float SampleSpacingXY)
{
	kernel_Localization << <dimGrid, dimBlock >> >(ParamIn, ParamNext, Convergence, FirstDev, SecondDev,
		Nfit, N_int, FitBoxCenter, Lambda, SampleSpacingXY);
}

extern "C" void kernel_getdev_wrapper(dim3 dimGrid, dim3 dimBlock, float *Data, float *gainR, float *PSFs, float *dPSFx, float *dPSFy, float *dPSFz, float *I, float *Bg, int Nfit, int PSFSize,
	float *FirstDev, float *SecondDev)
{
	kernel_getdev << <dimGrid, dimBlock >> >(Data, gainR, PSFs, dPSFx, dPSFy, dPSFz, I, Bg, Nfit, PSFSize, FirstDev, SecondDev);
}


extern "C" void kernel_calcPSFPixel_wrapper(dim3 dimGrid, dim3 dimBlock, float *SamplePSF, float *SamplePSFx, float *SamplePSFy, float *SamplePSFz, float *SamplePSFxy, float *SamplePSFxz, float *SamplePSFyz, float *SamplePSFxyz,
	int SizeX, int SizeY, float *PSFs, float *dPSFx, float *dPSFy, float *dPSFz, float * X, float *Y, float *Z,
	float SampleSpacingXY, float SampleSpacingZ, float StartX, float StartY, float StartZ, int N_int, int PSFSizeOut, int NPSFs)
{

	kernel_calcPSFPixel << <dimGrid, dimBlock >> >(SamplePSF, SamplePSFx, SamplePSFy, SamplePSFz, SamplePSFxy, SamplePSFxz, SamplePSFyz, SamplePSFxyz,
		SizeX, SizeY, PSFs, dPSFx, dPSFy, dPSFz, X, Y, Z,
		SampleSpacingXY, SampleSpacingZ, StartX, StartY, StartZ, N_int, PSFSizeOut, NPSFs);

}

extern "C" void kernel_calcPSFI_wrapper(dim3 dimGrid, dim3 dimBlock, float *psf, float *psfI, float *I, float *bg, int Nfit, int PSFsize)
{
	kernel_calcPSFI << <dimGrid, dimBlock >> >(psf, psfI, I, bg, Nfit, PSFsize);
}

extern "C" void kernel_Err_wrapper(dim3 dimGrid, dim3 dimBlock, float *PSF, float *Data, float *gainR, float *Error, int Nfit, int PSFSize)
{
	kernel_Err << <dimGrid, dimBlock >> >(PSF, Data, gainR, Error, Nfit, PSFSize);
}