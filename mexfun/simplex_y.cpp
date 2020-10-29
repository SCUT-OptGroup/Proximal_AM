#include "mex.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	double *Z, *b;

	double *theta;

	int m, nn;

	int k = 1;

	double sum1 = 0;

	mxAssert(nlhs == 1 && nrhs == 2, "Error: number of variables");

	m = mxGetM(prhs[0]);

	nn = mxGetN(prhs[0]);

	plhs[0] = mxCreateDoubleMatrix(1, nn, mxREAL);

	theta = (double*)mxGetData(plhs[0]);

	Z = (double*)mxGetData(prhs[0]);

	b = (double*)mxGetData(prhs[1]);


	for (int i = 1; i <= nn; i++)
	{

        sum1 = sum1 + Z[(i - 1)*m + k - 1];
        
		while ((Z[(i - 1)*m + k - 1] - (sum1 - b[i - 1]) / k > 0) && (k <= m - 1))
		{

			k = k + 1;
            
            sum1 = sum1 + Z[(i - 1)*m + k - 1];

		}

		if (Z[(i - 1)*m + k - 1] - (sum1 - b[i - 1]) / k <= 0)
		{

			sum1 = sum1 - Z[(i - 1)*m + k - 1];
            
            k = k - 1;

		}

		theta[i - 1] = (sum1 - b[i - 1]) / k;

		k = 1;

		sum1 = 0;

	}

}