// version: 1.0.
// @author: Jeppe Druedahl og Thomas Høgholm Jørgensen, 2016.

#include<cmath>
#include"mex.h"
#include"matrix.h"
#include<omp.h>
#include"base_funcs.c"
#include"mesh_interp.cpp"

#define MAX(X,Y) ((X)>(Y)?(X):(Y))
#define MIN(X,Y) ((X)<(Y)?(X):(Y))
#define MAXTHREADS MIN(16, omp_get_max_threads()-1)


/////////////
// GATEWAY //
/////////////

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    ////////////////////
    // 1. load inputs //
    ////////////////////

    int iin = 0; // first input

    // a. par: grids
    par_struct* par = new par_struct;

    par->ndim = (int) mxGetScalar(mxGetField(prhs[iin],0,"ndim"));

    par->grid_m = (double*) mxGetPr(mxGetField(prhs[iin],0,"grid_m"));
    par->grid_n = (double*) mxGetPr(mxGetField(prhs[iin],0,"grid_n"));
    par->grid_k = (double*) mxGetPr(mxGetField(prhs[iin],0,"grid_k"));
    par->Nm     = (int) mxGetScalar(mxGetField(prhs[iin],0,"Nm"));
    par->Nn     = (int) mxGetScalar(mxGetField(prhs[iin],0,"Nn"));
    par->Nk     = (int) mxGetScalar(mxGetField(prhs[iin],0,"Nk"));

    par->grid_a_pd = (double*) mxGetPr(mxGetField(prhs[iin],0,"grid_a_pd"));
    par->grid_b_pd = (double*) mxGetPr(mxGetField(prhs[iin],0,"grid_b_pd"));
    par->grid_q_pd = (double*) mxGetPr(mxGetField(prhs[iin],0,"grid_q_pd"));
    par->Na_pd     = (int)     mxGetScalar(mxGetField(prhs[iin],0,"Na_pd"));
    par->Nb_pd     = (int)     mxGetScalar(mxGetField(prhs[iin],0,"Nb_pd"));
    par->Nq_pd     = (int)     mxGetScalar(mxGetField(prhs[iin],0,"Nq_pd"));

    par->N_guess_vfi = (int) mxGetScalar(mxGetField(prhs[iin],0,"N_guess_vfi"));

    // b. par: parameters
    par->beta    = (double)  mxGetScalar(mxGetField(prhs[iin],0,"beta"));
    par->rho     = (double)  mxGetScalar(mxGetField(prhs[iin],0,"rho"));
    par->alpha   = (double)  mxGetScalar(mxGetField(prhs[iin],0,"alpha"));
    par->varphi  = (double)  mxGetScalar(mxGetField(prhs[iin],0,"varphi"));
    par->gamma   = (double)  mxGetScalar(mxGetField(prhs[iin],0,"gamma"));
    par->chi     = (double)  mxGetScalar(mxGetField(prhs[iin],0,"chi"));
    par->rk      = (double)  mxGetScalar(mxGetField(prhs[iin],0,"rk"));
    par->delta   = (double)  mxGetScalar(mxGetField(prhs[iin],0,"delta"));
    par->lmax    = (double)  mxGetScalar(mxGetField(prhs[iin],0,"lmax"));

    // c. par: settings
    par->max_threads     = (int)     mxGetScalar(mxGetField(prhs[iin],0,"max_threads"));

        // next intput
        iin++;

    // d. par: continuation value
    double* w = (double*) mxGetPr(mxGetField(prhs[iin],0,"w_values")); iin++;

    ////////////////////////
    // 2. allocate output //
    ////////////////////////

    int iout = 0; // first output

    // a. output
    int* dims = new int[par->ndim];
    dims[0]  = par->Nm;
    dims[1]  = par->Nn;
    if(par->ndim == 3){
        dims[2]  = par->Nk;
    }

    plhs[iout] = mxCreateNumericArray(par->ndim,dims,mxDOUBLE_CLASS,mxREAL);
    double* c  = mxGetPr(plhs[iout]); iout++;

    plhs[iout] = mxCreateNumericArray(par->ndim,dims,mxDOUBLE_CLASS,mxREAL);
    double* d  = mxGetPr(plhs[iout]); iout++;

    plhs[iout] = mxCreateNumericArray(par->ndim,dims,mxDOUBLE_CLASS,mxREAL);
    double* v  = mxGetPr(plhs[iout]); iout++;


    #pragma omp parallel num_threads(MAXTHREADS)
    {

        double** grid_vectors = new double*[par->ndim];
        grid_vectors[0] = par->grid_a_pd;
        grid_vectors[1] = par->grid_b_pd;

        int* grid_lengths = new int[par->ndim];
        grid_lengths[0] = par->Na_pd;
        grid_lengths[1] = par->Nb_pd;

        mesh::settings_struct* w_interp_func = new mesh::settings_struct;
        mesh::create(w_interp_func, par->ndim, grid_lengths, grid_vectors, w);

        double* pd = new double[par->ndim];

    // loop over states
    #pragma omp for
    for(int i_n =0; i_n < par->Nn; i_n++){
    for(int i_m =0; i_m < par->Nm; i_m++){

        double curr_c, curr_d;
        double curr_v = -mxGetInf();;

        // a. states
        double mi = par->grid_m[i_m];
        double ni = par->grid_n[i_n];

        // b. loop over guesses
        double min_c  = 0.001*mi;
        double max_c  = mi;
        double min_d  = 0.0;
        double max_d  = max_c;

        for(int i_cguess = 0; i_cguess < par->N_guess_vfi+1; i_cguess++){
        for(int i_dguess = 0; i_dguess < par->N_guess_vfi; i_dguess++){

            // i. c and d guess
            double c_guess, d_guess;
            if(i_cguess==par->N_guess_vfi){ // add points on the constraint
                double delta_c = (double)i_dguess/(double)(par->N_guess_vfi-1);
                c_guess = (1.0-delta_c)*min_c + delta_c*max_c;
                d_guess = mi + - c_guess;
            }  else {
                double delta_c = (double)i_cguess/(double)(par->N_guess_vfi-1);
                double delta_d = (double)i_dguess/(double)(par->N_guess_vfi-1);
                c_guess = (1.0-delta_c)*min_c + delta_c*max_c;
                d_guess = (1.0-delta_d)*min_d + delta_d*max_d;
            }

            // ii. check contraints
            if(c_guess + d_guess <= mi){

                // a. post-decision states
                double a_guess  = mi - c_guess - d_guess;
                double b_guess  = ni + d_guess + par->chi*log(1.0+d_guess);

                // b. continuation-value-of-choice
                pd[0] = a_guess;
                pd[1] = b_guess;
                double w_guess =  mesh::interp(w_interp_func,pd);

                // c. value-of-choice
                double v_guess = u(c_guess, 0, par->rho, par->alpha, par->varphi, par->gamma) + par->beta*w_guess;

                // d. update if highest value-of-choice
                if(v_guess > curr_v){
                    curr_v = v_guess;
                    curr_c = c_guess;
                    curr_d = d_guess;
                }

            } // allowed guess

        } // d
        } // c

        // c. fine-tune
        double delta_c = max_c*1.0/(double)(par->N_guess_vfi-1);
        min_c   = MAX(curr_c - delta_c,min_c);
        max_c   = MIN(curr_c + delta_c,max_c);

        double delta_d = max_d*1.0/(double)(par->N_guess_vfi-1);
        min_d   = MAX(curr_d - delta_d,min_d);
        max_d   = MIN(curr_d + delta_d,max_d);

        for(int i_cguess = 0; i_cguess < par->N_guess_vfi/4; i_cguess++){
        for(int i_dguess = 0; i_dguess < par->N_guess_vfi/4; i_dguess++){

            // i. c and d guess
            double delta_c = (double)i_cguess/(double)(par->N_guess_vfi/4-1);
            double delta_d = (double)i_dguess/(double)(par->N_guess_vfi/4-1);
            double c_guess = (1.0-delta_c)*min_c + delta_c*max_c;
            double d_guess = (1.0-delta_d)*min_d + delta_d*max_d;

            // ii. check contraints
            if(c_guess + d_guess <= mi){

                // a. post-decision states
                double a_guess  = mi - c_guess - d_guess;
                double b_guess  = ni + d_guess + par->chi*log(1.0+d_guess);

                // b. continuation-value-of-choice
                pd[0] = a_guess;
                pd[1] = b_guess;
                double w_guess =  mesh::interp(w_interp_func,pd);

                // c. value-of-choice
                double v_guess = u(c_guess, 0, par->rho, par->alpha, par->varphi, par->gamma) + par->beta*w_guess;

                // d. update if highest value-of-choice
                if(v_guess > curr_v){
                    curr_v = v_guess;
                    curr_c = c_guess;
                    curr_d = d_guess;
                }

            } // allowed guess

        } // d
        } // c


        // d. store results in the output matrix
        int igrid = index_func(i_m,i_n,0,par->Nm,par->Nn,1);
        v[igrid] = curr_v;
        c[igrid] = curr_c;
        d[igrid] = curr_d;

    } // i_n
    } // i_m

        delete[] pd;
        delete[] grid_lengths;
        delete[] grid_vectors;
        mesh::destroy(w_interp_func);


    } // end of parallel

} // GATEWAY
