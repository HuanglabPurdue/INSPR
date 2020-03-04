#include <windows.h>
#pragma comment(lib, "kernel32.lib")

//#include <iostream>
#include <mex.h>
#include "Model_Localization.hpp"


void CalNextParam(Model_Localization &obj, float *x0_next, float *x0, float *Converg, int Nfit);


void PopulArray(float *x, float *y, float *z, float *I, float *bg, float *err, int Nfit, float *x0);

// This function is mandatory.  It takes the place of 'main' in c.
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	/*!
 *  \brief Entry point in the code for Matlab.  Equivalent to main().
 *  \param nlhs number of left hand mxArrays to return
 *  \param plhs array of pointers to the output mxArrays
 *  \param nrhs number of input mxArrays
 *  \param prhs array of pointers to the input mxArrays.
 */
	//declare all vars
	// input data and paramters 
	float *channel1 = 0;
	float *gainR = 0;
	float *Coords1 = 0;
	float *xIn = 0;
	float *PSFsi = 0, *PSFxi = 0, *PSFyi = 0, *PSFzi = 0, *PSFxyi = 0, *PSFxzi = 0, *PSFyzi = 0, *PSFxyzi = 0; //spline parameters
	int  FitBoxSize = 16;
	int Nfit;
	float lambda;			// damping parameter in Marquardt-Levenberg algorithm
	float SampleSpacingXY, SampleSpacingZ;
	float StartX;			// Start of X dimention
	float StartY;			// Start of Y dimension
	float StartZ;			// Start of Z dimension

	const mwSize *datasize1 = 0, *datasize2 = 0, *datasize3 = 0;
	const mwSize *PSFsi_size = 0;
	mwSize outsize[1];
	mwSize outsize1[1];
	mwSize outsize2[1];
	mwSize outsize3[1];
	int iterateN1;
	float *result = 0;
	float *convergence = 0;
	// Local variables

	int Start, i, k;
	int SampledPSF_XYPixels, SampledPSF_ZPixels;
	float *x0 = 0;
	float *x0_next = 0;
	float *x0_Var = 0;
	int psfsize;
	int x0Size = NP;	//???
	int OTFpsize = 4;	//???
	int resultSize;
	float sse = 0;
	float LLR = 0;

	size_t NdimData1;
	//check for required inputs, correct types, and dimensions
	//1D vectors still return 2D
	datasize1 = mxGetDimensions(prhs[0]);
	datasize2 = mxGetDimensions(prhs[11]);

//	mexPrintf("datasize: %d %d\n", datasize1, datasize2);


	NdimData1 = mxGetNumberOfDimensions(prhs[0]);
	//retrieve all inputs	
	channel1 = (float *)mxGetData(prhs[0]);		//data; psfsize x Nfit x channel
	Coords1 = (float *)mxGetData(prhs[1]);
	PSFsi = (float*)mxGetData(prhs[2]);		//sample PSF; template

	SampleSpacingXY = (float)mxGetScalar(prhs[3]);		//XY pixel size in sample psf
	SampleSpacingZ = (float)mxGetScalar(prhs[4]);		//Z pixel size  in sample psf
	StartX = (float)mxGetScalar(prhs[5]);
	StartY = (float)mxGetScalar(prhs[6]);
	StartZ = (float)mxGetScalar(prhs[7]);

	iterateN1 = (int)mxGetScalar(prhs[8]);
	Nfit = (int)mxGetScalar(prhs[9]);		//Nfit
	lambda = (float)mxGetScalar(prhs[10]);
	xIn = (float*)mxGetData(prhs[11]);

	PSFxi = (float*)mxGetData(prhs[12]);	//spline parameters
	PSFyi = (float*)mxGetData(prhs[13]);
	PSFzi = (float*)mxGetData(prhs[14]);
	PSFxyi = (float*)mxGetData(prhs[15]);
	PSFxzi = (float*)mxGetData(prhs[16]);
	PSFyzi = (float*)mxGetData(prhs[17]);
	PSFxyzi = (float*)mxGetData(prhs[18]);

	gainR = (float*)mxGetData(prhs[19]);	//sCMOS noise



	//Caculated values
	const mwSize *SampledSize = mxGetDimensions(prhs[2]);
	SampledPSF_XYPixels = SampledSize[0];
	SampledPSF_ZPixels = SampledSize[2] / NCH;		//divide the number of channle; edited by FX


//	mexPrintf("test: %d %d\n", SampledSize[0], SampledSize[2]);

	//allocate memory for variables other than input

	x0 = (float *)mxCalloc(NPL*Nfit, sizeof(float));
	x0_next = (float *)mxCalloc(NPL*Nfit, sizeof(float));

	//validate input values(this section better not be blank!)
	//check dimensions
	if (NdimData1 != 2)
		mexErrMsgTxt("Array of channel1 must be a N x 1 vector.");
	//check sizes

	if (datasize1[1] != 1)
		mexErrMsgTxt("Array of channel1 must be a N x 1 vector.");
	if (datasize2[1] != 1)
		mexErrMsgTxt("Array of x0 must be a N x 1 vector.");
	if (FitBoxSize*FitBoxSize != datasize1[0] / NCH / Nfit)
		mexErrMsgTxt("boxsize must be square root of image size.");

	//check that input data are float

	if (!mxIsSingle(prhs[0]))
		mexErrMsgTxt("data must be type single");
	if (!mxIsSingle(prhs[1]))
		mexErrMsgTxt("Coords must be type single");
	if (!mxIsSingle(prhs[2]))
		mexErrMsgTxt("sample PSF must be type single");
	if (!mxIsSingle(prhs[11]))
		mexErrMsgTxt("x0 must be type single");

	for (i = 0; i < 7; i++)
	{
		if (!mxIsSingle(prhs[12 + i]))
			mexErrMsgTxt("sample PSF must be type single");
		PSFsi_size = mxGetDimensions(prhs[12 + i]);
		if (PSFsi_size[0] != SampledSize[0] || PSFsi_size[1] != SampledSize[1] || PSFsi_size[2] != SampledSize[2])
			mexErrMsgTxt("all sample PSF size must be the same");
	}

	// create output
	int errsize = 2;

	psfsize = FitBoxSize*FitBoxSize;
	outsize[0] = NPL*Nfit;
	outsize2[0] = errsize*Nfit;
	outsize3[0] = NCH*psfsize*Nfit;		//one channel, edited by FX
	plhs[0] = mxCreateNumericArray(1, outsize, mxSINGLE_CLASS, mxREAL);		//P_cuda
	plhs[1] = mxCreateNumericArray(1, outsize, mxSINGLE_CLASS, mxREAL);		//CG
	plhs[2] = mxCreateNumericArray(1, outsize, mxSINGLE_CLASS, mxREAL);		//crlb		
	plhs[3] = mxCreateNumericArray(1, outsize2, mxSINGLE_CLASS, mxREAL);	//err
	plhs[4] = mxCreateNumericArray(1, outsize3, mxSINGLE_CLASS, mxREAL);	//PSF_cuda
	result = (float *)mxGetData(plhs[0]);
	convergence = (float *)mxGetData(plhs[1]);
	//do stuff

	Start = (FitBoxSize % 2 == 0 ? FitBoxSize / 2 : (FitBoxSize - 1) / 2);
	memcpy(x0, xIn, NPL*Nfit*sizeof(float));
	//create PSF object and set properties 
	cudaDeviceReset();
	//create mLocalization object and set properties
	Model_Localization mFit;
	mFit.SampledPSF = PSFsi;
	mFit.SamplePSFx = PSFxi;	//spline parameters
	mFit.SamplePSFy = PSFyi;
	mFit.SamplePSFz = PSFzi;
	mFit.SamplePSFxy = PSFxyi;
	mFit.SamplePSFxz = PSFxzi;
	mFit.SamplePSFyz = PSFyzi;
	mFit.SamplePSFxyz = PSFxyzi;
	mFit.SampledPSF_XYPixels = SampledPSF_XYPixels;
	mFit.SampledPSF_ZPixels = SampledPSF_ZPixels;
	mFit.SampleSpacingXY = SampleSpacingXY;
	mFit.SampleSpacingZ = SampleSpacingZ;
	mFit.StartX = StartX;
	mFit.StartY = StartY;
	mFit.StartZ = StartZ;
	mFit.FitBoxCenter = Start;
	mFit.Data = channel1;
	mFit.Nfit = Nfit;
	mFit.Lambda = lambda;

	mFit.GainR = gainR;		//sCMOS noise


	mFit.ParamNext = (float*)malloc(NPL*Nfit*sizeof(float));
	mFit.ParamVar = (float *)mxGetData(plhs[2]);
	mFit.Convergence = (float*)malloc(NPL*Nfit*sizeof(float));
	mFit.PSFs = (float *)mxGetData(plhs[4]);
	mFit.Error = (float *)mxGetData(plhs[3]);	
	mFit.setup();
	//  iteration
	for (i = 0; i < iterateN1; i++)
	{
	//	std::cerr << "iteration: " << i << std::endl;

		CalNextParam(mFit, x0_next, x0, convergence, Nfit);

		//for (j = 0; j < x0Size*Nfit; j++) { x0[j] = x0_next[j]; }
		memcpy(x0, x0_next, NPL*Nfit*sizeof(float));
	}

	//write x0 in fitting result
	for (k = 0; k < NPL*Nfit; k++) result[k] = x0[k];
	

	// PSF Biplane, use x0,	 edited by FX
	float *xT = 0;
	float *yT = 0;
	float *zT = 0;
	float *IT = 0;
	float *bgT = 0;
	float *errT = 0;

	xT = (float*)malloc(Nfit*sizeof(float));
	yT = (float*)malloc(Nfit*sizeof(float));
	zT = (float*)malloc(Nfit*sizeof(float));
	IT = (float*)malloc(NCH * Nfit*sizeof(float));
	bgT = (float*)malloc(NCH * Nfit*sizeof(float));
	errT = (float*)malloc(2 * Nfit*sizeof(float));
	

	PopulArray(xT, yT, zT, IT, bgT, errT, Nfit, x0);

	//mFit.NPSFs = Nfit;
	mFit.sendParams(xT, yT, zT, IT, bgT, x0, errT);
	mFit.calPSF();
	mFit.calErr();
	mFit.calCRLB();
	
	free(xT);
	free(yT);
	free(zT);
	free(IT);
	free(bgT);
	free(mFit.ParamNext);
	free(mFit.Convergence);
	//free(mFit.ParamVar);

	//clean up anything not using an mxMalloc
	mxFree(x0);
	mxFree(x0_next);
	x0 = x0_next = 0;

	// We're done with the GPU now, so free the constant arrays
	mFit.cleanUp();
	cudaDeviceReset();
	return;
}

void CalNextParam(Model_Localization &obj, float *x0_next, float *x0, float *Converg, int Nfit)
{
	float *x = 0;
	float *y = 0;
	float *z = 0;
	float *I = 0;
	float *bg = 0;
	float *err = 0;

	x = (float*)mxMalloc(Nfit*sizeof(float));
	y = (float*)mxMalloc(Nfit*sizeof(float));
	z = (float*)mxMalloc(Nfit*sizeof(float));
	I = (float*)mxMalloc(NCH*Nfit*sizeof(float));
	bg = (float*)mxMalloc(NCH*Nfit*sizeof(float));
	err = (float*)mxMalloc(2 * Nfit*sizeof(float));

	PopulArray(x, y, z, I, bg, err, Nfit, x0);

	/* model generation and localization */

	obj.sendParams(x, y, z, I, bg, x0, err);
	obj.mLocalization();
	memcpy(x0_next, obj.ParamNext, NPL*Nfit*sizeof(float));
	memcpy(Converg, obj.Convergence, NPL*Nfit*sizeof(float));

	/* clean up */
	mxFree(x);
	mxFree(y);
	mxFree(z);
	mxFree(I);
	mxFree(bg);
	x = y = z = I = 0;
	bg = 0;

}


void PopulArray(float *x, float *y, float *z, float *I, float *bg, float *err, int Nfit, float *x0)
{
	int s, j;

	for (j = 0; j < Nfit; j++)
	{
		x[j] = x0[j*NPL];
		y[j] = x0[j*NPL + 1];
		z[j] = x0[j*NPL + 2];
		err[j * 2] = 0;
		err[j * 2 + 1] = 0;
		for (s = 0; s < NCH; s++)
		{

			I[s*Nfit + j] = x0[j*NPL + s + 3];	//??? need to change
			bg[s*Nfit + j] = x0[j*NPL + s + 3 + NCH];

		}
	}

}


