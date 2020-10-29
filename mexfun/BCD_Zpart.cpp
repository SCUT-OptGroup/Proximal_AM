#include "mex.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    double *Z, *ylambda, *n, *sig_alpha, *alpha, *beta, *WZ;
    
    double *temp_Zk;
    
    int m, nn, N;
    
    int sum1 = 0;
    
	mxAssert(nlhs == 1 && nrhs == 7, "Error: number of variables");

	m = mxGetM(prhs[0]);

	nn = mxGetN(prhs[0]);
    
    N = mxGetN(prhs[1]);
    
    plhs[0] = mxCreateDoubleMatrix(m, nn, mxREAL);

	temp_Zk = (double*)mxGetData(plhs[0]);

	Z = (double*)mxGetData(prhs[0]);

	ylambda = (double*)mxGetData(prhs[1]);

	n = (double*)mxGetData(prhs[2]);

	sig_alpha = (double*)mxGetData(prhs[3]);

	alpha = (double*)mxGetData(prhs[4]);

	beta = (double*)mxGetData(prhs[5]);

	WZ = (double*)mxGetData(prhs[6]);
    
    double *sum3 = new double[m*N];

    for (int s = 1; s <= m*N; s++)

    {
                
        sum3[s - 1] = 0;
                
    }
    
	for (int s = 1; s <= N; s++)
	{
        
        for (int j = 1; j <= n[s - 1]; j++)
		{

			for (int i = 1; i <= m; i++)
			{
                
				sum3[(s - 1)*m + i - 1] = sum3[(s - 1)*m + i - 1] + Z[(sum1 + j - 1)*m + i - 1];

			}

		}
        
		for (int j = 1; j <= n[s - 1]; j++)
		{

			for (int i = 1; i <= m; i++)
			{
                
				temp_Zk[(sum1 + j - 1)*m + i - 1] =  (1 / (sig_alpha[s - 1] + alpha[0])) * (sig_alpha[s - 1]*Z[(sum1 + j - 1)*m + i - 1] - beta[0]* sum3[(s - 1)*m + i - 1] + WZ[(sum1 + j - 1)*m + i - 1] + ylambda[(s - 1)*m + i - 1]);

			}

		}

        sum1 = sum1 + n[s-1];
            
	}
    
    delete []sum3;
    
}