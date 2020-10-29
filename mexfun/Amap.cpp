#include "mex.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	double *nn, *n, *X;

	double *Y;

	int m, N;

	int sum1 = 0;

	mxAssert(nlhs == 1 && nrhs == 3, "Error: number of variables");

	m = mxGetM(prhs[2]);

    N = mxGetN(prhs[2]);

	nn = (double*)mxGetData(prhs[0]);

	n = (double*)mxGetData(prhs[1]);

	X = (double*)mxGetData(prhs[2]);

	plhs[0] = mxCreateDoubleMatrix(m, nn[0], mxREAL);

	Y = (double*)mxGetData(plhs[0]);

	for (int s = 1; s <= m*N; s++)

	{

		Y[s - 1] = 0;

	}


	for (int s = 1; s <= N; s++)
	{

		for (int j = 1; j <= n[s - 1]; j++)
		{

			for (int i = 1; i <= m; i++)
			{

				Y[(sum1 + j - 1)*m + i - 1] = X[(s - 1)*m + i - 1];

			}

		}

		sum1 = sum1 + n[s - 1];

	}

}