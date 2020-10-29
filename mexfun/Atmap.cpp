#include "mex.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	double *N, *n, *X;

	double *Y;

	int m;

	int sum1 = 0;

	mxAssert(nlhs == 1 && nrhs == 3, "Error: number of variables");

	m = mxGetM(prhs[2]);

	N = (double*)mxGetData(prhs[0]);

	n = (double*)mxGetData(prhs[1]);

	X = (double*)mxGetData(prhs[2]);

	plhs[0] = mxCreateDoubleMatrix(m, N[0], mxREAL);

	Y = (double*)mxGetData(plhs[0]);

	for (int s = 1; s <= m*N[0]; s++)

	{

		Y[s - 1] = 0;

	}


	for (int s = 1; s <= N[0]; s++)
	{

		for (int j = 1; j <= n[s - 1]; j++)
		{

			for (int i = 1; i <= m; i++)
			{

				Y[(s - 1)*m + i - 1] = Y[(s - 1)*m + i - 1] + X[(sum1 + j - 1)*m + i - 1];

			}

		}

		sum1 = sum1 + n[s - 1];

	}

}