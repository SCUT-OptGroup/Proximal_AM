#include "mex.h"
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#define max(x,y) ( x>y?x:y )

int Compare(const void *a, const void *b)
{
    
    double da = *(double *)a; double db = *(double *)b;
    
    return (da < db) ? 1 : -1;
    
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	double *Z, *b; double *s_Z;

	int m, nn, i, j; int k = 1; double sum1 = 0;

	mxAssert(nlhs == 1 && nrhs == 2, "Error: number of variables");

	m = mxGetM(prhs[0]);

	nn = mxGetN(prhs[0]);

	plhs[0] = mxCreateDoubleMatrix(m, nn, mxREAL);
    
    s_Z = (double*)mxGetData(plhs[0]);

	Z = (double*)mxGetData(prhs[0]);

	b = (double*)mxGetData(prhs[1]);

    double** column = new double*[nn]; double* theta = new double[nn];

    #pragma omp parallel for   
	for (i = 1; i <= nn; i++)
	{
        
        column[i - 1] = new double[m];
        
        for (j = 1; j <= m; j++)
        {
            
            column[i - 1][j - 1] = Z[(i - 1)*m + j - 1];
        
        }
        
        qsort(column[i - 1],m,sizeof(column[i - 1][0]),Compare);
        
        for (j = 1; j <= m; j++)
        {
            
            s_Z[(i - 1)*m + j - 1] = column[i - 1][j - 1];
        
        }
        
        delete []column[i - 1];

    }
    
    delete []column;
    
    #pragma omp parallel for   
	for (int i = 1; i <= nn; i++)
	{
        
        sum1 = sum1 + s_Z[(i - 1)*m + k - 1];
        
		while ((s_Z[(i - 1)*m + k - 1] - (sum1 - b[i - 1]) / k > 0) && (k <= m - 1))
		{

			k = k + 1;
            
            sum1 = sum1 + s_Z[(i - 1)*m + k - 1];

		}

		if (s_Z[(i - 1)*m + k - 1] - (sum1 - b[i - 1]) / k <= 0)
		{

			sum1 = sum1 - s_Z[(i - 1)*m + k - 1];
            
            k = k - 1;

		}

		theta[i - 1] = (sum1 - b[i - 1]) / k;

		k = 1;

		sum1 = 0;

	}
    
      
	for (int i = 1; i <= nn; i++)
    {
        if (b[i - 1] == 0)  
        {
            for (j = 1; j <= m; j++)
            {
            
                s_Z[(i - 1)*m + j - 1] = 0;
        
            }
            
        }
        else
        {
            for (j = 1; j <= m; j++)
            {
            
                s_Z[(i - 1)*m + j - 1] = max(0,Z[(i - 1)*m + j - 1] - theta[i - 1]);
        
            }
            
        }
        
    }
    
    delete[] theta;
    
}
    
    