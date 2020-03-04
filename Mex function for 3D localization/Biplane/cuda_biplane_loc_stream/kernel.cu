
#include "cuda_runtime.h"
#include "definitions.h"
#include "kernel.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


__global__ void kernel_calcPSFPixel(float *SamplePSF, float *SamplePSFx, float *SamplePSFy, float *SamplePSFz, float *SamplePSFxy, float *SamplePSFxz, float *SamplePSFyz, float *SamplePSFxyz,
	int SizeX, int SizeY, float *PSFs, float *dPSFx, float *dPSFy, float *dPSFz, float * X, float *Y, float *Z,
	float SampleSpacingXY, float SampleSpacingZ, float StartX, float StartY, float StartZ, int N_int, int PSFSizeOut, int NPSFs, float *d_tMat, float *d_offset_int, int flag)

{
	__shared__ float PSF_Row_Samples[64], dPSFx_Row_Samples[64], dPSFy_Row_Samples[64], dPSFz_Row_Samples[64]; // dPSFx2_Row_Samples[64], dPSFy2_Row_Samples[64], dPSFz2_Row_Samples[64];	// store the interpolated pixels in here
	__shared__ float PSF_Row_Sum[16], dPSFx_Row_Sum[16], dPSFy_Row_Sum[16], dPSFz_Row_Sum[16];// dPSFx2_Row_Sum[16], dPSFy2_Row_Sum[16], dPSFz2_Row_Sum[16];		// This is the sum across rows and columns of interpolated pixels

	__shared__ float theta[5];				// X, Y, Z, I, Bg

	int idX = threadIdx.x;
	int PixelNumber = threadIdx.y;
	float F2d[2], dF2d_y[2], dF2d_x[2], dF2d_y2[2], dF2d_x2[2],
		dF2d_z[2], dF2d_zy[2], dF2d_zx[2], dF2d_zy2[2], dF2d_zx2[2], X1z[2];
	float dF[4];

	float a, b, c, d;	// coefficient for linear interpolation
	//Get XYZ, change to affine model
	//need change to affine position
	if ((threadIdx.x == 0) && (threadIdx.y == 0)){
		if (flag == 0) {
			theta[0] = X[blockIdx.y] * SampleSpacingXY * N_int;
			theta[1] = Y[blockIdx.y] * SampleSpacingXY * N_int;
		}
		else {

			float tmp1 = X[blockIdx.y] * SampleSpacingXY * N_int;
			float tmp2 = Y[blockIdx.y] * SampleSpacingXY * N_int;
			theta[0] = d_tMat[0] * tmp1 + d_tMat[1] * tmp2 + (d_tMat[2] + d_offset_int[blockIdx.y]) * SampleSpacingXY * N_int;
			theta[1] = d_tMat[3] * tmp1 + d_tMat[4] * tmp2 + (d_tMat[5] + d_offset_int[NPSFs + blockIdx.y]) * SampleSpacingXY * N_int;
	
		}

		theta[2] = Z[blockIdx.y];
	}

	__syncthreads();

	float Y_thread = (N_int*PixelNumber + idX + 0.5)*SampleSpacingXY - theta[1];

	int YBaseIndex = floor((Y_thread - StartX) / SampleSpacingXY);

	float SampleSpacingXYInv = 1 / SampleSpacingXY;
	float SampleSpacingZInv = 1 / SampleSpacingZ;
	int idZ = round((theta[2] - StartZ) * SampleSpacingZInv);


	//intialize PSF_Row_Sum counter
	if (threadIdx.x == 0)
	{
		PSF_Row_Sum[threadIdx.y] = 0;
		dPSFx_Row_Sum[threadIdx.y] = 0;
		dPSFy_Row_Sum[threadIdx.y] = 0;
		dPSFz_Row_Sum[threadIdx.y] = 0;

	}
	for (int ii = 0; ii < N_int; ii++) //go right in row
	{
		float X_thread = (blockIdx.x*N_int + ii + 0.5)*SampleSpacingXY - theta[0];

		//for interpolation we need the four surrounding points

		int XBaseIndex = floor((X_thread - StartY) * SampleSpacingXYInv);

		//using the follwing notation:
		//X1    X2
		//   o
		//X3	X4

		//These are values of the sampled points.

		//need change theta 0 and theta 1 with a certain offset  (X_thread, Y_thread)
		for (int nn = 0; nn < 2; nn++)
		{
			//These are locations of the sampled points
			float X1x = XBaseIndex*SampleSpacingXY + StartX;
			float X1y = YBaseIndex*SampleSpacingXY + StartY;
			X1z[nn] = (idZ + nn)*SampleSpacingZ + StartZ;

			int tmp2 = SizeY*SizeX*(idZ + nn) + SizeY*XBaseIndex + YBaseIndex;
			//Bicubic interpolation
			// first part
			interpXY(SamplePSF, SamplePSFx, SamplePSFy, SamplePSFxy, tmp2, SizeX, SampleSpacingXY, X1x, X1y,
				X_thread, Y_thread, &F2d[nn], &dF2d_y[nn], &dF2d_y2[nn], &dF2d_x[nn], &dF2d_x2[nn]);

			// second part
			interpXY(SamplePSFz, SamplePSFxz, SamplePSFyz, SamplePSFxyz, tmp2, SizeX, SampleSpacingXY, X1x, X1y,
				X_thread, Y_thread, &dF2d_z[nn], &dF2d_zy[nn], &dF2d_zy2[nn], &dF2d_zx[nn], &dF2d_zx2[nn]);
		}

		gencoeff(F2d[0], F2d[1], dF2d_z[0], dF2d_z[1], SampleSpacingZ, X1z[0], X1z[1], &a, &b, &c, &d);
		PSF_Row_Samples[N_int*PixelNumber + idX] = evalspline(a, b, c, d, theta[2]);
		dev_spline(a, b, c, theta[2], &dF[0]);
		dPSFz_Row_Samples[N_int*PixelNumber + idX] = dF[0];

		gencoeff(dF2d_x[0], dF2d_x[1], dF2d_zx[0], dF2d_zx[1], SampleSpacingZ, X1z[0], X1z[1], &a, &b, &c, &d);
		dPSFx_Row_Samples[N_int*PixelNumber + idX] = evalspline(a, b, c, d, theta[2]);

		gencoeff(dF2d_y[0], dF2d_y[1], dF2d_zy[0], dF2d_zy[1], SampleSpacingZ, X1z[0], X1z[1], &a, &b, &c, &d);
		dPSFy_Row_Samples[N_int*PixelNumber + idX] = evalspline(a, b, c, d, theta[2]);


		__syncthreads();
		//now sum over the row
		if (threadIdx.x == 0)
		for (int jj = 0; jj < N_int; jj++)
		{
			PSF_Row_Sum[threadIdx.y] += PSF_Row_Samples[N_int*PixelNumber + jj];
			dPSFx_Row_Sum[threadIdx.y] += dPSFx_Row_Samples[N_int*PixelNumber + jj];
			dPSFy_Row_Sum[threadIdx.y] += dPSFy_Row_Samples[N_int*PixelNumber + jj];
			dPSFz_Row_Sum[threadIdx.y] += dPSFz_Row_Samples[N_int*PixelNumber + jj];

		}
	}

	//now return value for each pixel

	__syncthreads();
	if (threadIdx.x == 0)
	{
		PSFs[PSFSizeOut*PSFSizeOut*blockIdx.y + PSFSizeOut*blockIdx.x + threadIdx.y] = PSF_Row_Sum[threadIdx.y];
		dPSFx[PSFSizeOut*PSFSizeOut*blockIdx.y + PSFSizeOut*blockIdx.x + threadIdx.y] = dPSFx_Row_Sum[threadIdx.y];
		dPSFy[PSFSizeOut*PSFSizeOut*blockIdx.y + PSFSizeOut*blockIdx.x + threadIdx.y] = dPSFy_Row_Sum[threadIdx.y];
		dPSFz[PSFSizeOut*PSFSizeOut*blockIdx.y + PSFSizeOut*blockIdx.x + threadIdx.y] = dPSFz_Row_Sum[threadIdx.y];

	}


}


__global__ void kernel_calcPSFI(float *psf, float *psfI, float *I, float *bg, int Nfit, int PSFsize)
{
	const int tx = threadIdx.x;
	const int bx = blockIdx.x;
	const int BlockSize = blockDim.x;

	//Prevent read/write past end of array
	int j = BlockSize*bx + tx;
	if ((bx*BlockSize + tx) >= Nfit) return;

	for (int i = 0; i < PSFsize; i++)
	{
		psfI[j*PSFsize + i] = psf[j*PSFsize + i] * I[j] + bg[j];

	}
}



__device__ void interpXY(float *SamplePSF, float *SamplePSFx, float *SamplePSFy, float *SamplePSFxy, int tmp2, int SizeX, float SampleSpacingXY, float X1x, float X1y,
	float X_thread, float Y_thread, float *F2d, float *dF2d_y, float *dF2d_y2, float *dF2d_x, float *dF2d_x2)
{
	float F[4], dF[4];
	float dX1X2[2], dX3X4[2], d2X1X2_y[2], d2X3X4_y[2];
	float a, b, c, d;	// coefficient for linear interpolation

	getsampleVal(SamplePSF, tmp2, SizeX, &F[0]);
	getsampleVal(SamplePSFx, tmp2, SizeX, &dF[0]);

	gencoeff(F[0], F[1], dF[0], dF[1], SampleSpacingXY, X1x, X1x + SampleSpacingXY, &a, &b, &c, &d);
	float X1X2 = evalspline(a, b, c, d, X_thread);
	dev_spline(a, b, c, X_thread, &dX1X2[0]);

	gencoeff(F[2], F[3], dF[2], dF[3], SampleSpacingXY, X1x, X1x + SampleSpacingXY, &a, &b, &c, &d);
	float X3X4 = evalspline(a, b, c, d, X_thread);
	dev_spline(a, b, c, X_thread, &dX3X4[0]);

	getsampleVal(SamplePSFy, tmp2, SizeX, &F[0]);
	getsampleVal(SamplePSFxy, tmp2, SizeX, &dF[0]);

	gencoeff(F[0], F[1], dF[0], dF[1], SampleSpacingXY, X1x, X1x + SampleSpacingXY, &a, &b, &c, &d);
	float dX1X2_y = evalspline(a, b, c, d, X_thread);
	dev_spline(a, b, c, X_thread, &d2X1X2_y[0]);

	gencoeff(F[2], F[3], dF[2], dF[3], SampleSpacingXY, X1x, X1x + SampleSpacingXY, &a, &b, &c, &d);
	float dX3X4_y = evalspline(a, b, c, d, X_thread);
	dev_spline(a, b, c, X_thread, &d2X3X4_y[0]);

	gencoeff(X1X2, X3X4, dX1X2_y, dX3X4_y, SampleSpacingXY, X1y, X1y + SampleSpacingXY, &a, &b, &c, &d);
	F2d[0] = evalspline(a, b, c, d, Y_thread);
	dev_spline(a, b, c, Y_thread, &dF[0]);
	dF2d_y[0] = dF[0];
	dF2d_y2[0] = dF[1];

	gencoeff(dX1X2[0], dX3X4[0], d2X1X2_y[0], d2X3X4_y[0], SampleSpacingXY, X1y, X1y + SampleSpacingXY, &a, &b, &c, &d);
	dF2d_x[0] = evalspline(a, b, c, d, Y_thread);
	gencoeff(dX1X2[1], dX3X4[1], d2X1X2_y[1], d2X3X4_y[1], SampleSpacingXY, X1y, X1y + SampleSpacingXY, &a, &b, &c, &d);
	dF2d_x2[0] = evalspline(a, b, c, d, Y_thread);
}

__device__ void gencoeff(float f1, float f2, float df1, float df2, float dx, float x1, float x2, float *a, float *b, float *c, float *d)
{
	a[0] = (df1 + df2) / dx / dx + 2 * (f1 - f2) / dx / dx / dx;
	b[0] = (df2 - df1) / 2 / dx - 3 * (x1 + x2)*a[0] / 2;
	c[0] = df1 - 3 * a[0] * x1*x1 - 2 * b[0] * x1;
	d[0] = f1 - a[0] * x1*x1*x1 - b[0] * x1*x1 - c[0] * x1;
}

__device__ float evalspline(float a, float b, float c, float d, float x)
{
	float f;
	f = x*x*x*a + x*x*b + x*c + d;
	return f;
}

__device__ void dev_spline(float a, float b, float c, float x, float *df)
{
	df[0] = 3 * a*x*x + 2 * b*x + c;
	df[1] = 6 * a*x + 2 * b;
}
__device__ void getsampleVal(float *SamplePSF, int ind, int SizeX, float *F)
{
	F[0] = SamplePSF[ind];
	F[1] = SamplePSF[ind + SizeX];
	F[2] = SamplePSF[ind + 1];
	F[3] = SamplePSF[ind + SizeX + 1];
}