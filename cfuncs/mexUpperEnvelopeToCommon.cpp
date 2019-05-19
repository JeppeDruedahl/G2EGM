// version: 1.0.
// @author: Jeppe Druedahl og Thomas Høgholm Jørgensen, 2016.

#include<cmath>
#include"mex.h"
#include"matrix.h"
#include<omp.h>
#include"base_funcs.c"
#include"mesh_interp.cpp"
#include"HighResTimer_class.hpp"

#define MAX(X,Y) ((X)>(Y)?(X):(Y))
#define MIN(X,Y) ((X)<(Y)?(X):(Y))

#define MAXTHREADS MIN(16, omp_get_max_threads()-1)
#define TIMER 0

typedef struct {

    double *m, *n, *k, *c, *d, *l; // states and choices
    bool *valid;                   // indicator
    int Na, Nb, Nq, num;           // dimension sizes and segment number

} seg_struct;

/////////////
// GATEWAY //
/////////////

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    int im_min, im_max, in_min, in_max, ik_min, ik_max;

        #if TIMER >= 1
        HighResTimer timer, timer_interp, timer_all;
        timer_all.StartTimer();
        double time_prep_gen    = 0.0;
        double time_prep_w      = 0.0;
        double time_interp      = 0.0;
        double time_holes       = 0.0;
        double time_copy        = 0.0;
        double time_interp_w    = 0.0;
        double time_interp_c    = 0.0;
        double time_interp_v    = 0.0;
        double time_interp_out  = 0.0;
        #endif // TIMER

    ////////////////////
    // 1. load inputs //
    ////////////////////

    int iin = 0; // first input

    // a. par: grids
    par_struct* par = new par_struct;

    par->ndim   = (int) mxGetScalar(mxGetField(prhs[iin],0,"ndim"));

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
    par->egm_extrap_add  = (int)     mxGetScalar(mxGetField(prhs[iin],0,"egm_extrap_add"));
    par->egm_extrap_w    = (double)  mxGetScalar(mxGetField(prhs[iin],0,"egm_extrap_w"));

        // next intput
        iin++;

    // d. par: continuation value
    double* w = (double*) mxGetPr(mxGetField(prhs[iin],0,"w_values")); iin++;

    // e. segment
    seg_struct* seg = new seg_struct;

    seg->m     = (double*) mxGetPr(mxGetField(prhs[iin],0,"m"));
    seg->n     = (double*) mxGetPr(mxGetField(prhs[iin],0,"n"));
    seg->k     = (double*) mxGetPr(mxGetField(prhs[iin],0,"k"));
    seg->c     = (double*) mxGetPr(mxGetField(prhs[iin],0,"c"));
    seg->d     = (double*) mxGetPr(mxGetField(prhs[iin],0,"d"));
    seg->l     = (double*) mxGetPr(mxGetField(prhs[iin],0,"l"));

    seg->valid = (bool*) mxGetData(mxGetField(prhs[iin],0,"valid"));

    seg->Na  = (int) mxGetScalar(mxGetField(prhs[iin],0,"Na"));
    seg->Nb  = (int) mxGetScalar(mxGetField(prhs[iin],0,"Nb"));
    seg->Nq  = (int) mxGetScalar(mxGetField(prhs[iin],0,"Nq"));
    seg->num = (int) mxGetScalar(mxGetField(prhs[iin],0,"num"));

        // next input
        iin++;

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

    plhs[iout]      = mxCreateNumericArray(par->ndim,dims,mxDOUBLE_CLASS,mxREAL);
    double* c_grid  = mxGetPr(plhs[iout]); iout++;

    plhs[iout]      = mxCreateNumericArray(par->ndim,dims,mxDOUBLE_CLASS,mxREAL);
    double* d_grid  = mxGetPr(plhs[iout]); iout++;

    double* l_grid;
    if(par->ndim == 3){
        plhs[iout]      = mxCreateNumericArray(par->ndim,dims,mxDOUBLE_CLASS,mxREAL);
        l_grid  = mxGetPr(plhs[iout]); iout++;
    }

    plhs[iout]      = mxCreateNumericArray(par->ndim,dims,mxDOUBLE_CLASS,mxREAL);
    double* v_grid  = mxGetPr(plhs[iout]); iout++;

    plhs[iout]      = mxCreateNumericArray(par->ndim,dims,mxLOGICAL_CLASS,mxREAL);
    bool* hole_grid = (bool*) mxGetData(plhs[iout]); iout++;

    // b. temporary container
    int noutput = par->ndim+1+1;
    int numel = par->Nm*par->Nn*par->Nk*noutput;
    int iv = 0;
    int ic = 1;
    int id = 2;
    int il, ihole;
    if(par->ndim == 3){
        il = 3;
        ihole = 4;
    } else {
        ihole = 3;
    }

    double* output_global  = new double[numel];
    double** output_locals = new double*[MIN(MAXTHREADS,par->max_threads)];

    for(int i = 0; i < MIN(MAXTHREADS,par->max_threads); i++){
        output_locals[i] = new double[numel];
    }

        for(int in = 0; in < par->Nn; in++){
        for(int ik = 0; ik < par->Nk; ik++){
        for(int im = 0; im < par->Nm; im++){

            int ioutput = noutput*index_func(im,in,ik,par->Nm,par->Nn,par->Nk);;

            output_global[ioutput + iv]    = -mxGetInf();
            output_global[ioutput + ihole] = 1;

        } // m
        } // n
        } // k

    /////////////
    // 3. main //
    /////////////

    #pragma omp parallel num_threads(MIN(MAXTHREADS,par->max_threads))
    {

        double* output = output_locals[omp_get_thread_num()];

        // i. initialize
        for(int in = 0; in < par->Nn; in++){
        for(int ik = 0; ik < par->Nk; ik++){
        for(int im = 0; im < par->Nm; im++){

            int ioutput = noutput*index_func(im,in,ik,par->Nm,par->Nn,par->Nk);;

            output[ioutput + iv]    = -mxGetInf();
            output[ioutput + ihole] = 1;

        } // m
        } // n
        } // k

        // ii. interpolant
        double* pd = new double[par->ndim];

        double** grid_vectors = new double*[par->ndim];
        grid_vectors[0] = par->grid_a_pd;
        grid_vectors[1] = par->grid_b_pd;
        if(par->ndim == 3){
            grid_vectors[2] = par->grid_q_pd;
        }

        int* grid_lengths = new int[par->ndim];
        grid_lengths[0] = par->Na_pd;
        grid_lengths[1] = par->Nb_pd;
        if(par->ndim == 3){
            grid_lengths[2] = par->Nq_pd;
        }

        mesh::settings_struct* w_interp_func = new mesh::settings_struct;
        mesh::create(w_interp_func, par->ndim, grid_lengths, grid_vectors, w);


    // main loop
    #pragma omp for schedule(dynamic)
    for(int ib = 0; ib < seg->Nb; ib++){
    for(int iq = 0; iq < seg->Nq; iq++){
    for(int ia = 0; ia < seg->Na; ia++){

        int tri_num = 2;
        if(par->ndim == 3){
            tri_num = 4;
        }

    for(int tri = 0; tri < tri_num; tri++){

            #if TIMER >= 1
            timer.StartTimer();
            #endif // TIMER

        // a. simplex in (a,b,q)-space (or similar with constrained choices)
        int i1, i2, i3, i4;

            // fixed point
            i1 = index_func(ia,ib,iq,seg->Na,seg->Nb,seg->Nq);

        if(par->ndim == 3){ // 3d

            if(tri == 0){

                if(ia == seg->Na-1){continue;}
                if(ib == seg->Nb-1){continue;}
                if(iq == seg->Nq-1){continue;}

                i2 = index_func(ia+1,ib,iq,seg->Na,seg->Nb,seg->Nq);
                i3 = index_func(ia,ib+1,iq,seg->Na,seg->Nb,seg->Nq);
                i4 = index_func(ia,ib,iq+1,seg->Na,seg->Nb,seg->Nq);

            } else if(tri == 1){

                if(ia == seg->Na-1){continue;}
                if(ib == seg->Nb-1){continue;}
                if(iq == 0){continue;}

                i2 = index_func(ia+1,ib,iq,seg->Na,seg->Nb,seg->Nq);
                i3 = index_func(ia,ib+1,iq,seg->Na,seg->Nb,seg->Nq);
                i4 = index_func(ia,ib,iq-1,seg->Na,seg->Nb,seg->Nq);

            } else  if(tri == 2){

                if(ia == 0){continue;}
                if(ib == 0){continue;}
                if(iq == seg->Nq-1){continue;}

                i2 = index_func(ia-1,ib,iq,seg->Na,seg->Nb,seg->Nq);
                i3 = index_func(ia,ib-1,iq,seg->Na,seg->Nb,seg->Nq);
                i4 = index_func(ia,ib,iq+1,seg->Na,seg->Nb,seg->Nq);

            } else if(tri == 3){

                if(ia == 0){continue;}
                if(ib == 0){continue;}
                if(iq == 0){continue;}

                i2 = index_func(ia-1,ib,iq,seg->Na,seg->Nb,seg->Nq);
                i3 = index_func(ia,ib-1,iq,seg->Na,seg->Nb,seg->Nq);
                i4 = index_func(ia,ib,iq-1,seg->Na,seg->Nb,seg->Nq);

            } else if(tri == 4){

                if(ia == seg->Na-1){continue;}
                if(ib == 0){continue;}
                if(iq == seg->Nq-1){continue;}

                i2 = index_func(ia+1,ib,iq,seg->Na,seg->Nb,seg->Nq);
                i3 = index_func(ia,ib-1,iq,seg->Na,seg->Nb,seg->Nq);
                i4 = index_func(ia,ib,iq+1,seg->Na,seg->Nb,seg->Nq);

            } else if(tri == 5){

                if(ia == seg->Na-1){continue;}
                if(ib == 0){continue;}
                if(iq == 0){continue;}

                i2 = index_func(ia+1,ib,iq,seg->Na,seg->Nb,seg->Nq);
                i3 = index_func(ia,ib-1,iq,seg->Na,seg->Nb,seg->Nq);
                i4 = index_func(ia,ib,iq-1,seg->Na,seg->Nb,seg->Nq);

            } else  if(tri == 6){

                if(ia == 0){continue;}
                if(ib == seg->Nb-1){continue;}
                if(iq == seg->Nq-1){continue;}

                i2 = index_func(ia-1,ib,iq,seg->Na,seg->Nb,seg->Nq);
                i3 = index_func(ia,ib+1,iq,seg->Na,seg->Nb,seg->Nq);
                i4 = index_func(ia,ib,iq+1,seg->Na,seg->Nb,seg->Nq);

            } else if(tri == 7){

                if(ia == 0){continue;}
                if(ib == seg->Nb-1){continue;}
                if(iq == 0){continue;}

                i2 = index_func(ia-1,ib,iq,seg->Na,seg->Nb,seg->Nq);
                i3 = index_func(ia,ib+1,iq,seg->Na,seg->Nb,seg->Nq);
                i4 = index_func(ia,ib,iq-1,seg->Na,seg->Nb,seg->Nq);

            }

            // b. exit if one of the corners is not a valid choice
            if(seg->valid[i1] == false ||
               seg->valid[i2] == false ||
               seg->valid[i3] == false ||
               seg->valid[i4] == false){
                continue;
            }

        } else { // 2d

            // a. triangle in (a,b)-space (or similar with constrained choices)
            i1 = index_func(ia,ib,iq,seg->Na,seg->Nb,seg->Nq);
            if(ib == seg->Nb-1){continue;}
            i2 = index_func(ia,ib+1,iq,seg->Na,seg->Nb,seg->Nb);

            if(tri == 0){
                if(ia == 0){continue;}
                if(ib == seg->Nb-1){continue;}
                i3 = index_func(ia-1,ib+1,iq,seg->Na,seg->Nb,seg->Nq);
            } else {
                if(ia == seg->Na-1){continue;}
                i3 = index_func(ia+1,ib,iq,seg->Na,seg->Nb,seg->Nq);
            }

            // b. exit if one of the corners is not a valid choice
            if(seg->valid[i1] == false ||
               seg->valid[i2] == false ||
               seg->valid[i3] == false) {
                continue;
            }

        }

        // c. simplex in (m,n,k)-space
        double m1 = seg->m[i1];
        double m2 = seg->m[i2];
        double m3 = seg->m[i3];

        double n1 = seg->n[i1];
        double n2 = seg->n[i2];
        double n3 = seg->n[i3];

        double m4, n4, k1, k2, k3, k4;
        if(par->ndim == 3){

            m4 = seg->m[i4];
            n4 = seg->n[i4];
            k1 = seg->k[i1];
            k2 = seg->k[i2];
            k3 = seg->k[i3];
            k4 = seg->k[i4];
        }

        // d. boundary box values and indices
        double m_max = MAX(m1,MAX(m2,m3));
        double m_min = MIN(m1,MIN(m2,m3));
        double n_max = MAX(n1,MAX(n2,n3));
        double n_min = MIN(n1,MIN(n2,n3));

        double k_max, k_min;
        if(par->ndim == 3){
            k_max = MAX(k1,MAX(k2,k3));
            k_min = MIN(k1,MIN(k2,k3));
        }

        int im_low, in_low, ik_low;
        if(m_min < 0){
            im_low = 0;
        } else {
            im_low  = mesh::binary_search(0, par->Nm, &par->grid_m[0], m_min);
        }
        int im_high = mesh::binary_search(0, par->Nm, &par->grid_m[0], m_max) + 1;
        if(n_min < 0){
            in_low = 0;
        } else {
            in_low  = mesh::binary_search(0, par->Nn, &par->grid_n[0], n_min);
        }
        int in_high = mesh::binary_search(0, par->Nn, &par->grid_n[0], n_max) + 1;
        int ik_high;
        if(par->ndim == 3){
            if(k_min < 0){
                ik_low = 0;
            } else {
                ik_low  = mesh::binary_search(0, par->Nk, &par->grid_k[0], k_min);
            }
            ik_high = mesh::binary_search(0, par->Nk, &par->grid_k[0], k_max) + 1;
        }

            // correction to allow for more extrapolation
            im_low  = MAX(im_low-par->egm_extrap_add,0);
            im_high = MIN(im_high+par->egm_extrap_add,par->Nm);
            in_low  = MAX(in_low-par->egm_extrap_add,0);
            in_high = MIN(in_high+par->egm_extrap_add,par->Nn);
            if(par->ndim == 3){
                ik_low  = MAX(ik_low-par->egm_extrap_add,0);
                ik_high = MIN(ik_high+par->egm_extrap_add,par->Nk);
            }

        // e. prepare barycentric interpolation

            #if TIMER >= 1
            time_prep_gen += timer.StopTimer();
            timer.StartTimer();
            #endif // TIMER

        double denom, inv_detT, fac11, fac12, fac13, fac21, fac22, fac23, fac31, fac32, fac33;
        if(par->ndim == 3){

            double x_1 = m1-m4;
            double x_2 = m2-m4;
            double x_3 = m3-m4;
            double y_1 = n1-n4;
            double y_2 = n2-n4;
            double y_3 = n3-n4;
            double z_1 = k1-k4;
            double z_2 = k2-k4;
            double z_3 = k3-k4;

            double detT = -x_3*y_2*z_1+x_2*y_3*z_1+x_3*y_1*z_2-x_1*y_3*z_2-x_2*y_1*z_3+x_1*y_2*z_3;
            inv_detT = 1.0/detT;

            fac11 = (y_2*z_3-y_3*z_2);
            fac12 = (x_3*z_2-x_2*z_3);
            fac13 = (x_2*y_3-x_3*y_2);

            fac21 = (y_3*z_1-y_1*z_3);
            fac22 = (x_1*z_3-x_3*z_1);
            fac23 = (x_3*y_1-x_1*y_3);

            fac31 = (y_1*z_2-y_2*z_1);
            fac32 = (x_2*z_1-x_1*z_2);
            fac33 = (x_1*y_2-x_2*y_1);

        } else { // 2d

            denom = (n2-n3)*(m1-m3)+(m3-m2)*(n1-n3);

        }

            #if TIMER >= 1
            time_prep_w += timer.StopTimer();
            timer_interp.StartTimer();
            #endif // TIMER

        // f. loop through common grid nodes in interior of bounding box
        if(par->ndim == 2){
            ik_low = 0;
            ik_high = 1;
        }
        for(int ik = ik_low; ik < ik_high; ik++){
        for(int in = in_low; in < in_high; in++){
        for(int im = im_low; im < im_high; im++){

                #if TIMER >= 1
                timer.StartTimer();
                #endif // TIMER

            // i. common grid values
            double m_now = par->grid_m[im];
            double n_now = par->grid_n[in];

            double k_now;
            double w1, w2, w3, w4;
            if(par->ndim == 3){

                k_now = par->grid_k[ik];

                // ii. barycentric coordinates
                w1 = ( fac11*(m_now-m4) + fac12*(n_now-n4) + fac13*(k_now-k4) ) * inv_detT;
                w2 = ( fac21*(m_now-m4) + fac22*(n_now-n4) + fac23*(k_now-k4) ) * inv_detT;
                w3 = ( fac31*(m_now-m4) + fac32*(n_now-n4) + fac33*(k_now-k4) ) * inv_detT;
                w4 = 1.0 - w1 - w2 - w3;

                    #if TIMER >= 1
                    time_interp_w += timer.StopTimer();
                    #endif // TIMER

                // iii. exit if too much outside simplex
                if(w1 < par->egm_extrap_w ||
                   w2 < par->egm_extrap_w ||
                   w3 < par->egm_extrap_w ||
                   w4 < par->egm_extrap_w){
                        continue;
                }

            } else {

                // ii. barycentric coordinates
                w1 = ((n2-n3)*(m_now-m3) + (m3-m2)*(n_now-n3)) / denom;
                w2 = ((n3-n1)*(m_now-m3) + (m1-m3)*(n_now-n3)) / denom;
                w3 = 1 - w1 - w2;

                // iii. exit if too much outside simplex
                if(w1 < par->egm_extrap_w ||
                   w2 < par->egm_extrap_w ||
                   w3 < par->egm_extrap_w){
                        continue;
                }

            }

                #if TIMER >= 1
                timer.StartTimer();
                #endif // TIMER

            // iv. interpolate choices
            double c_interp, d_interp, l_interp, a_interp, b_interp, q_interp;
            if(par->ndim == 3){ // 3d
            if(seg->num == 1 || seg->num == 5){ // ucon and lcon, interpolate c, d, and maybe l

                c_interp = w1*seg->c[i1] + w2*seg->c[i2] +  w3*seg->c[i3] + w4*seg->c[i4];
                d_interp = w1*seg->d[i1] + w2*seg->d[i2] +  w3*seg->d[i3] + w4*seg->d[i4];
                if(seg->num == 1){
                    l_interp = w1*seg->l[i1] + w2*seg->l[i2] +  w3*seg->l[i3] + w4*seg->l[i4];
                } else {
                    l_interp = par->lmax;
                }
                a_interp  = m_now + par->rk*k_now*l_interp - c_interp - d_interp;
                b_interp  = n_now + d_interp + par->chi*log(1.0+d_interp);
                q_interp  = (1.0-par->delta)*k_now + l_interp;

            } else if(seg->num == 2 || seg->num == 6) { // dcon and dlcon, interpolate c and maybe l


                c_interp = w1*seg->c[i1] + w2*seg->c[i2] +  w3*seg->c[i3] + w4*seg->c[i4];
                d_interp = 0;
                if(seg->num == 2){
                    l_interp = w1*seg->l[i1] + w2*seg->l[i2] +  w3*seg->l[i3] + w4*seg->l[i4];
                } else {
                    l_interp = par->lmax;
                }

                a_interp  = m_now + par->rk*k_now*l_interp - c_interp - d_interp;
                b_interp  = n_now; // d_interp = 0
                q_interp  = (1.0-par->delta)*k_now + l_interp;

            } else if(seg->num == 3 || seg->num == 7) { // acon and lacon, interpolate d and maybe l

                a_interp = 0.0;
                d_interp = w1*seg->d[i1] + w2*seg->d[i2] +  w3*seg->d[i3] + w4*seg->d[i4];
                if(seg->num == 3){
                    l_interp = w1*seg->l[i1] + w2*seg->l[i2] +  w3*seg->l[i3] + w4*seg->l[i4];
                } else {
                    l_interp = par->lmax;
                }

                c_interp = m_now + par->rk*k_now*l_interp - a_interp - d_interp;
                b_interp = n_now + d_interp + par->chi*log(1.0+d_interp);
                q_interp  = (1.0-par->delta)*k_now + l_interp;

            } else if(seg->num == 4 || seg->num == 8) { // con or fullcon, interpolate l or nothing

                a_interp = 0.0;
                d_interp = 0.0;
                if(seg->num == 4){
                    l_interp = w1*seg->l[i1] + w2*seg->l[i2] +  w3*seg->l[i3] + w4*seg->l[i4];
                } else {
                    l_interp = par->lmax;
                }

                c_interp = m_now + par->rk*k_now*l_interp - a_interp - d_interp;
                b_interp = n_now + d_interp + par->chi*log(1.0+d_interp);
                q_interp  = (1.0-par->delta)*k_now + l_interp;

            }
            } else { // 2d

                l_interp = 0;
                q_interp = 0;

                if(seg->num == 1){ // ucon, interpolate c and d

                    c_interp = w1*seg->c[i1] + w2*seg->c[i2] + w3*seg->c[i3];
                    d_interp = w1*seg->d[i1] + w2*seg->d[i2] + w3*seg->d[i3];
                    a_interp  = m_now - c_interp - d_interp;
                    b_interp  = n_now + d_interp + par->chi*log(1.0+d_interp);

                } else if(seg->num == 2) { // dcon, interpolate c

                    c_interp = w1*seg->c[i1] + w2*seg->c[i2] + w3*seg->c[i3];
                    d_interp = 0.0;
                    a_interp  = m_now - c_interp - d_interp;
                    b_interp  = n_now; // d_interp = 0

                } else if(seg->num == 3) { // acon, interpolate d

                    a_interp = 0.0;
                    d_interp = w1*seg->d[i1] + w2*seg->d[i2] + w3*seg->d[i3];
                    c_interp = m_now - a_interp - d_interp;
                    b_interp = n_now + d_interp + par->chi*log(1.0+d_interp);

                }

            }
                #if TIMER >= 1
                time_interp_c += timer.StopTimer();
                #endif // TIMER

                // exit if choices are illegal
                if(c_interp <= 0.0 || d_interp < 0.0 || l_interp < 0.0){ continue; }
                if(a_interp < 0.0 || b_interp < 0.0 || q_interp < 0.0){ continue; }

            // v. value-of-choice

                #if TIMER >= 1
                timer.StartTimer();
                #endif // TIMER

            pd[0] = a_interp;
            pd[1] = b_interp;
            if(par->ndim == 3){
                pd[2] = q_interp;
            }

            double w_interp = mesh::interp(w_interp_func,pd);
            double v_interp = u(c_interp, l_interp, par->rho, par->alpha, par->varphi, par->gamma) + par->beta*w_interp;

                #if TIMER >= 1
                time_interp_v += timer.StopTimer();
                timer.StartTimer();
                #endif // TIMER

            // vi. update if max
            int ioutput = noutput*index_func(im,in,ik,par->Nm,par->Nn,par->Nk);

            if(v_interp > output[ioutput + iv]){

                output[ioutput + iv] = v_interp;

                output[ioutput + ic] = c_interp;
                output[ioutput + id] = d_interp;
                if(par->ndim == 3){
                    output[ioutput + il] = l_interp;
                }
                output[ioutput + ihole] = 0; // no hole


            } // max

                #if TIMER >= 1
                time_interp_out += timer.StopTimer();
                #endif // TIMER

        } // m
        } // n
        } // k

            #if TIMER >= 1
            time_interp += timer_interp.StopTimer();
            #endif // TIMER


    } // tri

    } // a
    } // b
    } // q

    ///////////////////////
    // 4. copy to global //
    ///////////////////////

    #pragma omp for schedule(dynamic)
    for(int in = 0; in < par->Nn; in++){
    for(int ik = 0; ik < par->Nk; ik++){
    for(int im = 0; im < par->Nm; im++){

        int igrid = index_func(im,in,ik,par->Nm,par->Nn,par->Nk);
        int ioutput = noutput*igrid;

        int imax = -1;
        double vmax;

        for(int i = 0; i < MIN(MAXTHREADS,par->max_threads); i++){
            if(output_locals[i][ioutput + ihole] == 0){
                if(imax == -1){
                    imax = i;
                    vmax = output_locals[i][ioutput + iv];
                } else if(output_locals[i][ioutput + iv] > vmax) {
                    imax = i;
                    vmax = output_locals[i][ioutput + iv];
                }
            }
        }

        if(imax > -1){

            output_global[ioutput + ihole] = 0;

            output_global[ioutput + iv] = output_locals[imax][ioutput + iv];
            output_global[ioutput + ic] = output_locals[imax][ioutput + ic];
            output_global[ioutput + id] = output_locals[imax][ioutput + id];

            if(par->ndim == 3){
                output_global[ioutput + il] = output_locals[imax][ioutput + il];
            }
        }

    } // n
    } // k
    } // m

    delete[] output_locals[omp_get_thread_num()];
    output = output_global;


    /////////////////////////////
    // 5. fill out empty holes //
    /////////////////////////////

        #if TIMER >= 1
        timer.StartTimer();
        #endif // TIMER

    #pragma omp single
    {

    // a. locate global bounding box
    in_min = 0;
    in_max = par->Nn-1;
    double min_n  = mxGetInf();
    double max_n  = -mxGetInf();

    im_min = 0;
    im_max = par->Nm-1;
    double min_m  = mxGetInf();
    double max_m  = -mxGetInf();

    ik_min = 0;
    ik_max = par->Nk-1;
    double min_k, max_k;
    if(par->ndim == 3){
        min_k = mxGetInf();
        max_k = -mxGetInf();
    } else {
        ik_min = 0;
        ik_max = 0;
    }

    for(int ik = 0; ik < par->Nk; ik++){
    for(int in = 0; in < par->Nn; in++){
    for(int im = 0; im < par->Nm; im++){

        int ioutput = noutput*index_func(im,in,ik,par->Nm,par->Nn,par->Nk);

        double m_now = par->grid_m[im];
        double n_now = par->grid_n[in];
        double k_now;
        if(par->ndim == 3){
            k_now = par->grid_k[ik];
        }

        if(m_now < min_m && (bool)output[ioutput + ihole] == 0){
            min_m  = m_now;
            im_min = im;
        }
        if(m_now > max_m && (bool)output[ioutput + ihole] == 0){
            max_m  = m_now;
            im_max = im;
        }

        if(n_now < min_n && (bool)output[ioutput + ihole] == 0){
            min_n  = n_now;
            in_min = in;
        }
        if(n_now > max_n && (bool)output[ioutput + ihole] == 0){
            max_n  = n_now;
            in_max = in;
        }

        if(par->ndim == 3){
        if(k_now < min_k && (bool)output[ioutput + ihole] == 0){
            min_k  = k_now;
            ik_min = ik;
        }
        if(k_now > max_k && (bool)output[ioutput + ihole] == 0){
            max_k  = k_now;
            ik_max = ik;
        }
        }

    } // k
    } // n
    } // m
    } // master
    #pragma omp barrier // wait until done

    // b. loop through m, n, k nodes to detect holes
    #pragma omp for schedule(dynamic)
    for(int in = in_min; in < MIN(in_max+1,par->Nn); in++){
    for(int ik = ik_min; ik < MIN(ik_max+1,par->Nk); ik++){
    for(int im = im_min; im < MIN(im_max+1,par->Nm); im++){

        int ioutput = noutput*index_func(im,in,ik,par->Nm,par->Nn,par->Nk);
        if(output[ioutput + ihole] == 0 ){ continue; }

        double m_now = par->grid_m[im];
        double n_now = par->grid_n[in];
        double k_now;
        int k_add = 0;
        int m_add = 2;
        int n_add = 2;
        if(par->ndim == 3){
            k_now = par->grid_k[ik];
            m_add = 3;
            n_add = 3;
            k_add = 3;
        }

        // loop over points close by
        for(int ik_close = MAX(0,ik-k_add); ik_close <= MIN(ik+k_add,par->Nk-1); ik_close++){
        for(int in_close = MAX(0,in-n_add); in_close <= MIN(in+n_add,par->Nn-1); in_close++){
        for(int im_close = MAX(0,im-m_add); im_close <= MIN(im+m_add,par->Nm-1); im_close++){

            int ioutput_close = noutput*index_func(im_close,in_close,ik_close,par->Nm,par->Nn,par->Nk);

                // not useable if hole itself
                if((bool)output[ioutput_close + ihole] == true){continue;}

            // i. copy choices from close point (nearest neighbor)
            double c_interp, d_interp, l_interp, a_interp, b_interp, q_interp;
            if(par->ndim == 3){ // 3d
            if(seg->num == 1 || seg->num == 5){ // ucon, interpolate c, d and l

                c_interp = output[ioutput_close + ic];
                d_interp = output[ioutput_close + id];
                l_interp = output[ioutput_close + il];

                a_interp  = m_now + par->rk*k_now*l_interp - c_interp - d_interp;
                b_interp  = n_now + d_interp + par->chi*log(1.0+d_interp);
                q_interp  = (1.0-par->delta)*k_now + l_interp;

            } else if(seg->num == 2 || seg->num == 6) { // dcon, interpolate c and l

                c_interp = output[ioutput_close + ic];
                d_interp = 0.0;
                l_interp = output[ioutput_close + il];

                a_interp  = m_now + par->rk*k_now*l_interp - c_interp - d_interp;
                b_interp  = n_now; // d_interp = 0
                q_interp  = (1.0-par->delta)*k_now + l_interp;

            } else if(seg->num == 3 || seg->num == 7) { // acon, interpolate d and l

                a_interp = 0.0;
                d_interp = output[ioutput_close + id];
                l_interp = output[ioutput_close + il];

                c_interp = m_now + par->rk*k_now*l_interp - a_interp - d_interp;
                b_interp = n_now + d_interp + par->chi*log(1.0+d_interp);
                q_interp  = (1.0-par->delta)*k_now + l_interp;

            } else if(seg->num == 4 || seg->num == 8) { // con, interpolate l

                a_interp = 0.0;
                d_interp = 0.0;
                l_interp = output[ioutput_close + il];

                c_interp = m_now + par->rk*k_now*l_interp - a_interp - d_interp;
                b_interp = n_now + d_interp + par->chi*log(1.0+d_interp);
                q_interp  = (1.0-par->delta)*k_now + l_interp;

            }
            } else {

                l_interp = 0;
                q_interp = 0;

                if(seg->num == 1){ // ucon, interpolate c and d

                    c_interp = output[ioutput_close + ic];
                    d_interp = output[ioutput_close + id];
                    a_interp  = m_now - c_interp - d_interp;
                    b_interp  = n_now + d_interp + par->chi*log(1.0+d_interp);

                } else if(seg->num == 2) { // dcon, interpolate c

                    c_interp = output[ioutput_close + ic];
                    d_interp = 0.0;
                    a_interp  = m_now - c_interp - d_interp;
                    b_interp  = n_now; // d_interp = 0

                } else if(seg->num == 3) { // acon, interpolate d

                    a_interp = 0.0;
                    d_interp = output[ioutput_close + id];
                    c_interp = m_now - a_interp - d_interp;
                    b_interp = n_now + d_interp + par->chi*log(1.0+d_interp);

                }
            }

               // exist if choices are illegal
               if(c_interp <= 0.0 || d_interp < 0.0 || l_interp < 0.0){ continue; }
               if(a_interp < 0.0 || b_interp < 0.0 || q_interp < 0.0){ continue; }

            // ii. value-of-choice
            pd[0] = a_interp;
            pd[1] = b_interp;
            if(par->ndim == 3){
                pd[2] = q_interp;
            }

            double w_interp = mesh::interp(w_interp_func,pd);
            double v_interp = u(c_interp, l_interp, par->rho, par->alpha, par->varphi, par->gamma) + par->beta*w_interp;

            // iii. save if highest value
            if(v_interp >  output[ioutput + iv]){

                output[ioutput + iv] = v_interp;

                output[ioutput + ic] = c_interp;
                output[ioutput + id] = d_interp;
                if(par->ndim == 3){
                    output[ioutput + il] = l_interp;
                }
            }

            } // m_close
            } // n_close
            } // k_close

    } // m
    } // k
    } // n

        #if TIMER >= 1
        time_holes += timer.StopTimer();
        timer.StartTimer();
        #endif // TIMER

    ///////////////////////
    // 6. copy to output //
    ///////////////////////

    #pragma omp for schedule(dynamic)
    for(int in = 0; in < par->Nn; in++){
    for(int ik = 0; ik < par->Nk; ik++){
    for(int im = 0; im < par->Nm; im++){

        int igrid = index_func(im,in,ik,par->Nm,par->Nn,par->Nk);
        int ioutput = noutput*igrid;

        v_grid[igrid] = output[ioutput + iv];

        if(mxIsFinite(v_grid[igrid]) == false){
             c_grid[igrid] = mxGetNaN();
             d_grid[igrid] = mxGetNaN();
             if(par->ndim == 3){
                l_grid[igrid] = mxGetNaN();
             }
        } else {
            c_grid[igrid] = output[ioutput + ic];
            d_grid[igrid] = output[ioutput + id];
            if(par->ndim == 3){
                l_grid[igrid] = output[ioutput + il];
            }
        }

        hole_grid[igrid] = (bool)output[ioutput + ihole];

    } // m
    } // k
    } // n

        #if TIMER >= 1
        time_copy += timer.StopTimer();
        timer.StartTimer();
        #endif // TIMER

        // clean up
        delete[] pd;
        delete[] grid_lengths;
        delete[] grid_vectors;
        mesh::destroy(w_interp_func);


    } // parallel: upper envelope

    // clean up
    delete seg;
    delete par;
    delete[] dims;
    delete[] output_locals;
    delete[] output_global;

        #if TIMER >= 1
        printf("    all:      %4.2f secs\n",timer_all.StopTimer());
        printf("     prep:    %4.2f secs\n",time_prep_gen);
        printf("     prep_w:  %4.2f secs\n",time_prep_w);
        printf("     interp:  %4.2f secs\n",time_interp);
        printf("      w:      %4.2f secs\n",time_interp_w);
        printf("      c:      %4.2f secs\n",time_interp_c);
        printf("      v:      %4.2f secs\n",time_interp_v);
        printf("      out:    %4.2f secs\n",time_interp_v);
        printf("     holes:   %4.2f secs\n",time_holes);
        printf("     copy:    %4.2f secs\n",time_copy);
        #endif // TIMER

} // GATEWAY
