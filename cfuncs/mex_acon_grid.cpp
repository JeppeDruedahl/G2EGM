// version: 1.0.
// @author: Jeppe Druedahl og Thomas Høgholm Jørgensen, 2016.

#include<cmath>
#include"mex.h"
#include"matrix.h"

void nonlinespace(double *x, double low, double high, int n, double phi)
{
    int i;

    x[0] = low;
    for(i = 1; i < n-1; i++){
        x[i] = x[i-1] + (high - x[i-1])/(pow((double)(n-i),phi));
    }
    x[n-1] = high; // numerically exact

    return;

} // nonlinespace

/////////////
// GATEWAY //
/////////////

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // input
    int iin = 0;
    double* c_min   = (double*) mxGetPr(prhs[iin]); iin++;
    double* c_max   = (double*) mxGetPr(prhs[iin]); iin++;

    double* grid_b_acon = (double*) mxGetPr(mxGetField(prhs[iin],0,"grid_b_acon")); // vector
    double* grid_q_acon = (double*) mxGetPr(mxGetField(prhs[iin],0,"grid_q_acon")); // vector
    int Nc_acon         = (int)     mxGetScalar(mxGetField(prhs[iin],0,"Nc_acon"));
    int Nb_acon         = (int)     mxGetScalar(mxGetField(prhs[iin],0,"Nb_acon"));
    int Nq_acon         = (int)     mxGetScalar(mxGetField(prhs[iin],0,"Nq_acon"));
    double phi_m        = (double)  mxGetScalar(mxGetField(prhs[iin],0,"phi_m"));

    // output
    int dims[3];
    dims[0]  = Nc_acon;
    dims[1]  = Nb_acon;
    dims[2]  = Nq_acon;

    int iout = 0;
    plhs[iout]      = mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);
    double* c_acon  = mxGetPr(plhs[iout]); iout++;
    plhs[iout]      = mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);
    double* b_acon  = mxGetPr(plhs[iout]); iout++;
    plhs[iout]      = mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);
    double* q_acon  = mxGetPr(plhs[iout]); iout++;

    // fill output
    double* vec = new double[Nc_acon];
    for(int ib = 0; ib < Nb_acon; ib++){
    for(int iq = 0; iq < Nq_acon; iq++){

        double c_min_now = c_min[iq*Nb_acon+ib];
        double c_max_now = c_max[iq*Nb_acon+ib];

        nonlinespace(vec, c_min_now, c_max_now, Nc_acon, 1.0);

        for(int ic = 0; ic < Nc_acon; ic++){
            c_acon[iq*Nb_acon*Nc_acon + ib*Nc_acon + ic] = vec[ic];
            b_acon[iq*Nb_acon*Nc_acon + ib*Nc_acon + ic] = grid_b_acon[ib];
            q_acon[iq*Nb_acon*Nc_acon + ib*Nc_acon + ic] = grid_q_acon[iq];
        }
    } // q
    } // b

    // free memory
    delete[] vec;

}
