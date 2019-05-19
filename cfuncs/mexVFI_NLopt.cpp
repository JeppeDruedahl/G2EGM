// version: 1.0.
// @author: Jeppe Druedahl og Thomas Høgholm Jørgensen, 2016.

#include<cmath>
#include"mex.h"
#include"matrix.h"
#include<omp.h>
#include"base_funcs.c"
#include"mesh_interp.cpp"
#include"nlopt-2.4.2-dll64\nlopt.h"

#define MAX(X,Y) ((X)>(Y)?(X):(Y))
#define MIN(X,Y) ((X)<(Y)?(X):(Y))
#define MAXTHREADS MIN(16, omp_get_max_threads()-1)

///////////////////
// SOLVER STRUCT //
///////////////////

typedef struct
{

    par_struct *par;                      // parameters
    double mi, ni, ki;                    // states
    mesh::settings_struct* w_interp_func; // interpolant
    double *pd;
    double *w;

} solver_struct;


////////////////////////
// OBJECTIVE FUNCTION //
////////////////////////

double v_guess(double c_guess, double d_guess, double l_guess, solver_struct *solver_data)
{

    par_struct* par = solver_data->par;

    double mi = solver_data->mi;
    double ni = solver_data->ni;
    double ki;
    if(par->ndim == 3){
        ki = solver_data->ki;
    }

    // 1. post-decision states
    double a_guess, b_guess, q_guess;
    if(par->ndim == 3){ // 3d

        a_guess  = mi + par->rk*ki*l_guess - c_guess - d_guess;
        b_guess  = ni + d_guess + par->chi*log(1.0+d_guess);
        q_guess  = (1.0-par->delta)*ki + l_guess;

        solver_data->pd[0] = a_guess;
        solver_data->pd[1] = b_guess;
        solver_data->pd[2] = q_guess;

    } else { // 2d

        a_guess  = mi - c_guess - d_guess;
        b_guess  = ni + d_guess + par->chi*log(1.0+d_guess);

        solver_data->pd[0] = a_guess;
        solver_data->pd[1] = b_guess;

    }
    double w_guess = mesh::interp(solver_data->w_interp_func,solver_data->pd);
    double v_guess = u(c_guess, l_guess, par->rho, par->alpha, par->varphi, par->gamma) + par->beta*w_guess;

    return -v_guess;

}

double objfunc(unsigned n, const double *x, double *grad, void *solver_data_in)
{

    solver_struct *solver_data = (solver_struct *) solver_data_in;
    double forward, eps=1e-5;

    // 1. choices
    double c_guess = x[0];
    double d_guess = x[1];
    double l_guess;
    if(solver_data->par->ndim == 3){
        l_guess = x[2];
    }

    // 2. value of choice
    double obj = v_guess(c_guess, d_guess, l_guess, solver_data);

    if (grad) {

        forward = v_guess(c_guess + eps, d_guess, l_guess, solver_data);
        grad[0] = (forward - obj)/eps;

        forward = v_guess(c_guess, d_guess + eps, l_guess, solver_data);
        grad[1] = (forward - obj)/eps;

        if(solver_data->par->ndim == 3){
            forward = v_guess(c_guess, d_guess, l_guess + eps, solver_data);
            grad[2] = (forward - obj)/eps;
        }
    }
    return obj;

}


////////////////////////////
// INEQUALITY CONSTRAINTS //
////////////////////////////

double ineq_constraint(unsigned n, const double *x, double *grad, void *solver_data_in)
{
    solver_struct *solver_data = (solver_struct *) solver_data_in;
    par_struct* par = solver_data->par;

    if (grad) {
        grad[0] = 1.0;
        grad[1] = 1.0;
        if(par->ndim == 3){
            grad[2] = -par->rk*solver_data->ki;
        }
    }
    if(par->ndim == 3){
        return (x[0] + x[1]) - (solver_data->mi + par->rk*solver_data->ki*x[2]); // positive -> violated
    } else {
        return x[0] + x[1] - solver_data->mi; // positive -> violated
    }
 }


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
    size_t* dims = new size_t[par->ndim];
    dims[0]  = par->Nm;
    dims[1]  = par->Nn;
    if(par->ndim == 3){
        dims[2]  = par->Nk;
    }

    plhs[iout] = mxCreateNumericArray(par->ndim,dims,mxDOUBLE_CLASS,mxREAL);
    double* c  = mxGetPr(plhs[iout]); iout++;

    plhs[iout] = mxCreateNumericArray(par->ndim,dims,mxDOUBLE_CLASS,mxREAL);
    double* d  = mxGetPr(plhs[iout]); iout++;

    double* l;
    if(par->ndim == 3){
        plhs[iout] = mxCreateNumericArray(par->ndim,dims,mxDOUBLE_CLASS,mxREAL);
        l  = mxGetPr(plhs[iout]); iout++;
    }

    plhs[iout] = mxCreateNumericArray(par->ndim,dims,mxDOUBLE_CLASS,mxREAL);
    double* v  = mxGetPr(plhs[iout]); iout++;


    /////////////
    // 3. main //
    /////////////

    #pragma omp parallel num_threads(MIN(MAXTHREADS,par->max_threads))
    {


    // a. setup solver data
    solver_struct *solver_data = new solver_struct;

        // parameters
        solver_data->par = par;

        // interpolant
        double** grid_vectors = new double*[par->ndim];
        grid_vectors[0] = par->grid_a_pd;
        grid_vectors[1] = par->grid_b_pd;
        if(par->ndim == 3){
            grid_vectors[2] = par->grid_q_pd;
        }

        int* grid_lengths    = new int[par->ndim];
        grid_lengths[0] = par->Na_pd;
        grid_lengths[1] = par->Nb_pd;
        if(par->ndim == 3){
            grid_lengths[2] = par->Nq_pd;
        }

        solver_data->w_interp_func = new mesh::settings_struct;
        mesh::create(solver_data->w_interp_func, par->ndim, grid_lengths, grid_vectors, w);

        solver_data->pd = new double[par->ndim];
        solver_data->w = w;

    // b. create optimization object;
    nlopt_opt opt;
    opt = nlopt_create(NLOPT_LD_MMA, par->ndim);

        double* x = new double[par->ndim];
        double* lb = new double[par->ndim];
        double* ub = new double[par->ndim];

        // settings
        nlopt_set_min_objective(opt, objfunc, solver_data);
        nlopt_set_xtol_rel(opt, 1e-6);
        nlopt_set_ftol_rel(opt, 1e-6);
        nlopt_set_maxeval(opt, 500);

        // bounds
        lb[0] = 0;
        lb[1] = 0;
        if(par->ndim == 3){
            lb[2] = 0;
        }
        nlopt_set_lower_bounds(opt, lb);

        if(par->lmax > 0){
            ub[0] = mxGetInf();
            ub[1] = mxGetInf();
            ub[2] = par->lmax;
            nlopt_set_upper_bounds(opt, ub);
        }

        // constraints
        nlopt_add_inequality_constraint(opt, ineq_constraint, solver_data, 1e-8);

    // c. loop over states
    #pragma omp for
    for(int i_m = 0; i_m < par->Nm; i_m++){
    for(int i_k = 0; i_k < par->Nk; i_k++){
    for(int i_n = 0; i_n < par->Nn; i_n++){

    int igrid = index_func(i_m,i_n,i_k,par->Nm,par->Nn,par->Nk);
    v[igrid]  = -mxGetInf();

    // multistart loop
    int ini_max;
    if(par->ndim == 3){
        if(i_n == 0){
            ini_max = 6;
        } else {
            ini_max = 7;
        }
    } else {
        if(i_n == 0){
            ini_max = 3;
        } else {
            ini_max = 4;
        }
    }
    for(int ini = 0; ini < ini_max; ini++){

        // i. states
        double mi = par->grid_m[i_m];
        double ni = par->grid_n[i_n];
        double ki;
        if(par->ndim == 3){
            ki = par->grid_k[i_k];
        }
            // update solver data
            solver_data->mi = mi;
            solver_data->ni = ni;
            solver_data->ki = ki;

            if(par->ndim == 2){
                ub[0] = mi;
                ub[1] = mi;
                nlopt_set_upper_bounds(opt, ub);
            }

        // ii. initial guess
        if(par->ndim ==3 ){ // 3d
        if(ini < 3){

            if(i_n == 0){
                x[2] = 0.1*par->lmax;
            } else {
                int igrid_lag = index_func(i_m,i_n-1,i_k,par->Nm,par->Nn,par->Nk);
                x[2] = 0.9*l[igrid_lag];
            }
            if(ini == 0){
                // inner conner
                x[0] = (mi+par->rk*ki*x[2])*1.0/5.0;
                x[1] = (mi+par->rk*ki*x[2])*1.0/5.0;
            } else if(ini == 1){
                // lower right corner
                x[0] = (mi+par->rk*ki*x[2])*4.0/5.0;
                x[1] = (mi+par->rk*ki*x[2])*1.0/6.0;
            } else if(ini == 2){
                // upper corner
                x[0] = (mi+par->rk*ki*x[2])*1.0/6.0;
                x[1] = (mi+par->rk*ki*x[2])*4.0/5.0;
            }

        } else if(ini < 6) {

            if(i_n == 0){
                x[2] = 0.99*par->lmax;
            } else {
                int igrid_lag = index_func(i_m,i_n-1,i_k,par->Nm,par->Nn,par->Nk);
                x[2] = 1.1*l[igrid_lag];
            }
            if(ini == 3){
                // inner conner
                x[0] = (mi+par->rk*ki*x[2])*1.0/5.0;
                x[1] = (mi+par->rk*ki*x[2])*1.0/5.0;
            } else if(ini == 4){
                // lower right corner
                x[0] = (mi+par->rk*ki*x[2])*4.0/5.0;
                x[1] = (mi+par->rk*ki*x[2])*1.0/6.0;
            } else if(ini == 5){
                // upper corner
                x[0] = (mi+par->rk*ki*x[2])*1.0/6.0;
                x[1] = (mi+par->rk*ki*x[2])*4.0/5.0;
            }

        } else {

            int igrid_lag = index_func(i_m,i_n-1,i_k,par->Nm,par->Nn,par->Nk);
            x[0] = c[igrid_lag];
            x[1] = d[igrid_lag];
            x[2] = l[igrid_lag];

        }
        } else { // 2d

            if(ini == 0){
                // inner conner
                x[0] = mi*1.0/5.0;
                x[1] = mi*1.0/5.0;
            } else if(ini == 1){
                // lower right corner
                x[0] = mi*4.0/5.0;
                x[1] = mi*1.0/6.0;
            } else if(ini == 2){
                // upper corner
                x[0] = mi*1.0/6.0;
                x[1] = mi*4.0/5.0;
            } else {
                int igrid_lag = index_func(i_m,i_n-1,i_k,par->Nm,par->Nn,par->Nk);
                x[0] = c[igrid_lag];
                x[1] = d[igrid_lag];
            }

        }

        // iii. run solver
        double minf;
        nlopt_optimize(opt, x, &minf);

        // iv. update solution
        if(-minf > v[igrid]){

            v[igrid] = -minf;

            c[igrid] = x[0];
            d[igrid] = x[1];

            if(par->ndim == 3){
                l[igrid] = x[2];
            }

        }

    } // ini;
    } // m
    } // n
    } // k

        nlopt_destroy(opt);
        delete[] grid_lengths;
        delete[] grid_vectors;
        mesh::destroy(solver_data->w_interp_func);
        delete[] x;
        delete[] lb;
        delete[] ub;
        delete[] solver_data->pd;
        delete solver_data;

    } // end of parallel

    // clean up
    delete par;
    delete[] dims;

} // GATEWAY
