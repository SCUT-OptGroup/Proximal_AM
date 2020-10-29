#include "mex.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    double *Z, *y, *lambda, *tau, *sigma, *beta, *alpha, *N, *n;

	double *dinf1;

	int m, nn;

	int sum1 = 0;

	mxAssert(nlhs == 1 && nrhs == 9, "Error: number of variables");

	m = mxGetM(prhs[0]);

	nn = mxGetN(prhs[0]);

	plhs[0] = mxCreateDoubleMatrix(m, nn, mxREAL);

	dinf1 = (double*)mxGetData(plhs[0]);

	Z = (double*)mxGetData(prhs[0]);

	y = (double*)mxGetData(prhs[1]);

	lambda = (double*)mxGetData(prhs[2]);

	tau = (double*)mxGetData(prhs[3]);

	sigma = (double*)mxGetData(prhs[4]);

	beta = (double*)mxGetData(prhs[5]);

	alpha = (double*)mxGetData(prhs[6]);

	N = (double*)mxGetData(prhs[7]);
    
    n = (double*)mxGetData(prhs[8]);
    
    double *sum3 = new double[m*N[0]];
    
    for (int s = 1; s <= m*N[0]; s++)

    {
                
        sum3[s - 1] = 0;
                
    }

	for (int s = 1; s <= N[0]; s++)
	{

		for (int j = 1; j <= n[s - 1]; j++)
		{

			for (int i = 1; i <= m; i++)
			{

				sum3[(s - 1)*m + i - 1] = sum3[(s - 1)*m + i - 1] + Z[(sum1 + j - 1)*m + i - 1];

			}

		}

		sum1 = sum1 + n[s - 1];

	}

	sum1 = 0;

	for (int s = 1; s <= N[0]; s++)
	{

		for (int j = 1; j <= n[s - 1]; j++)
		{

			for (int i = 1; i <= m; i++)
			{

				dinf1[(sum1 + j - 1)*m + i - 1] = sigma[s - 1]*Z[(sum1 + j - 1)*m + i - 1] + beta[0] * (y[i - 1] - sum3[(s - 1)*m + i - 1]) + ((1 - tau[0])/tau[0])*lambda[(s - 1)*m + i - 1];

			}

		}

		sum1 = sum1 + n[s - 1];

	}
    
    delete []sum3;
    
}