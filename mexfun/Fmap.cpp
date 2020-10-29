#include "mex.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	double *X, *a;

	double *W;

	int m, nn, d;

	double sum1 = 0, temp;

	mxAssert(nlhs == 1 && nrhs == 2, "Error: number of variables");

	d = mxGetM(prhs[0]);

	m = mxGetN(prhs[0]);

	nn = mxGetN(prhs[1]);

	plhs[0] = mxCreateDoubleMatrix(m, nn, mxREAL);

	W = (double*)mxGetData(plhs[0]);

	X = (double*)mxGetData(prhs[0]);

	a = (double*)mxGetData(prhs[1]);


	for (int j = 1; j <= nn; j++)
	{

		for (int i = 1; i <= m; i++)
		{

			for (int t = 1; t <= d; t++)
			{

				temp = X[(i - 1)*d + t - 1] - a[(j - 1)*d + t - 1];

				sum1 = sum1 + temp*temp;

			}

			W[(j - 1)*m + i - 1] = sum1;

			sum1 = 0;

		}

	}

}