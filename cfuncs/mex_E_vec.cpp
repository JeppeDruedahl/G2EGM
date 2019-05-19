// version: 1.0.
// @author: Jeppe Druedahl og Thomas Høgholm Jørgensen, 2016.

#include<cmath>
#include"mex.h"
#include"matrix.h"
#include<omp.h>
#include"base_funcs.c"
#include"mesh_interp.cpp"
#include"vectorclass\vectorclass.h"
#include"vectorclass\vectormath_exp.h"

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

    par->eta             = (double*) mxGetPr(mxGetField(prhs[iin],0,"eta"));
    par->w_eta           = (double*) mxGetPr(mxGetField(prhs[iin],0,"w_eta"));
    par->Neta            = (int)     mxGetScalar(mxGetField(prhs[iin],0,"Neta"));
    par->Nm_ret          = (int)     mxGetScalar(mxGetField(prhs[iin],0,"Nm_ret"));
    par->do_derivatives  = (int)     mxGetScalar(mxGetField(prhs[iin],0,"do_derivatives"));

    // b. par: parameters
    par->yret        = (double)  mxGetScalar(mxGetField(prhs[iin],0,"yret"));
    par->Ra          = (double)  mxGetScalar(mxGetField(prhs[iin],0,"Ra"));
    par->Rb          = (double)  mxGetScalar(mxGetField(prhs[iin],0,"Rb"));
    par->rk_retire   = (double)  mxGetScalar(mxGetField(prhs[iin],0,"rk_retire"));
    par->sigma       = (double)  mxGetScalar(mxGetField(prhs[iin],0,"sigma"));
    double inv_sigma = 1.0/par->sigma;

    // c. par: settings
    par->max_threads     = (int)     mxGetScalar(mxGetField(prhs[iin],0,"max_threads"));

        // next intput
        iin++;

    // sol: v and vprime
    double *v = (double*) mxGetPr(prhs[iin]); iin++;
    double *vm = (double*) mxGetPr(prhs[iin]); iin++;
    double *vn = (double*) mxGetPr(prhs[iin]); iin++;
    double *vk;
    if(par->ndim == 3){
        vk = (double*) mxGetPr(prhs[iin]); iin++;
    }

    double* grid_m_retire  = (double*) mxGetPr(mxGetField(prhs[iin],0,"m"));
    double* v_retire       = (double*) mxGetPr(mxGetField(prhs[iin],0,"v"));
    double *vm_retire, *vn_retire, *vk_retire;
    if(par->do_derivatives == 1){
        vm_retire      = (double*) mxGetPr(mxGetField(prhs[iin],0,"vm"));
        vn_retire      = (double*) mxGetPr(mxGetField(prhs[iin],0,"vn"));
        if(par->ndim == 3){
            vk_retire = (double*) mxGetPr(mxGetField(prhs[iin],0,"vk"));
        }
    }

    ////////////////////////
    // 2. allocate output //
    ////////////////////////

    int iout = 0; // first output

    // a. output
    size_t* dims = new size_t[par->ndim];
    dims[0]  = par->Na_pd;
    dims[1]  = par->Nb_pd;
    if(par->ndim == 3){
        dims[2]  = par->Nq_pd;
    }

    plhs[iout]       = mxCreateNumericArray(par->ndim,dims,mxDOUBLE_CLASS,mxREAL);
    double* w_grid   = mxGetPr(plhs[iout]); iout++;

    plhs[iout]       = mxCreateNumericArray(par->ndim,dims,mxDOUBLE_CLASS,mxREAL);
    double* wa_grid  = mxGetPr(plhs[iout]); iout++;

    plhs[iout]       = mxCreateNumericArray(par->ndim,dims,mxDOUBLE_CLASS,mxREAL);
    double* wb_grid  = mxGetPr(plhs[iout]); iout++;

    double* wq_grid;
    if(par->ndim == 3){
        plhs[iout]       = mxCreateNumericArray(par->ndim,dims,mxDOUBLE_CLASS,mxREAL);
        double* wq_grid  = mxGetPr(plhs[iout]); iout++;
    }

    /////////////
    // 3. main //
    /////////////

    printf("main\n");
    #pragma omp parallel num_threads(MIN(MAXTHREADS,par->max_threads))
    {

        // i. interpolation prep
        int dimf = 1;
        if(par->do_derivatives == 1){
            dimf = 1 + par->ndim;
        }
        double** yi = new double*[dimf];

        // ii. interpolant work
        double** v_vprime         = new double*[dimf];
        for(int i = 0; i < dimf; i++){
            v_vprime[i] = new double[par->Neta];
        }
        double* next_states = new double[par->ndim];

        double** grids_v_vprime   = new double*[dimf];
        double** grid_vectors = new double*[par->ndim];
        grid_vectors[0] = par->grid_m;
        grid_vectors[1] = par->grid_n;
        if(par->ndim == 3){
            grid_vectors[2] = par->grid_k;
        }

        int* grid_lengths = new int[par->ndim];
        grid_lengths[0] = par->Nm;
        grid_lengths[1] = par->Nn;
        if(par->ndim == 3){
            grid_lengths[2] = par->Nk;
        }

        mesh::settings_struct* v_vprime_interp_func = new mesh::settings_struct;
        if(par->do_derivatives == 0){
            grids_v_vprime[0] = v;
            mesh::create_multiple(v_vprime_interp_func, par->ndim, grid_lengths, grid_vectors, dimf, grids_v_vprime);
        } else if(par->ndim == 3){
            grids_v_vprime[0] = v;
            grids_v_vprime[1] = vm;
            grids_v_vprime[2] = vn;
            grids_v_vprime[3] = vk;
            mesh::create_multiple(v_vprime_interp_func, par->ndim, grid_lengths, grid_vectors, dimf, grids_v_vprime);
        } else {
            grids_v_vprime[0] = v;
            grids_v_vprime[1] = vm;
            grids_v_vprime[2] = vn;
            mesh::create_multiple(v_vprime_interp_func, par->ndim, grid_lengths, grid_vectors, dimf, grids_v_vprime);
        }

        // iii. interpolant retire
        double** v_vprime_retire      = new double*[dimf];
        for(int i = 0; i < dimf; i++){
            v_vprime_retire[i] = new double[par->Neta];
        }

        double** grid_v_vprime_retire = new double*[dimf];
        double** grid_vector_retire = new double*[1];
        grid_vector_retire[0] = grid_m_retire;

        int* grid_length_retire = new int[1];
        grid_length_retire[0] = par->Nm_ret;

        mesh::settings_struct* v_vprime_retire_interp_func = new mesh::settings_struct;
        if(par->do_derivatives == 0){
            grid_v_vprime_retire[0] = v_retire;
            mesh::create_multiple(v_vprime_retire_interp_func, 1, grid_length_retire, grid_vector_retire, dimf, grid_v_vprime_retire);
        } else if(par->ndim == 3){
            grid_v_vprime_retire[0] = v_retire;
            grid_v_vprime_retire[1] = vm_retire;
            grid_v_vprime_retire[2] = vn_retire;
            grid_v_vprime_retire[3] = vk_retire;
            mesh::create_multiple(v_vprime_retire_interp_func, 1, grid_length_retire, grid_vector_retire, dimf, grid_v_vprime_retire);
        } else {
            grid_v_vprime_retire[0] = v_retire;
            grid_v_vprime_retire[1] = vm_retire;
            grid_v_vprime_retire[2] = vn_retire;
            mesh::create_multiple(v_vprime_retire_interp_func, 1, grid_length_retire, grid_vector_retire, dimf, grid_v_vprime_retire);
        }


        // iv. position vectors
        int* pos_left = new int[par->ndim];
        int* pos_left_retire = new int[1];
        pos_left[1] = 0; // n_next is always increasing (b is first loop, never affected but eta)

        // v. eta vectors
        Vec4d eta_vec, w_eta_vec;
        w_eta_vec.load(par->w_eta);
        eta_vec.load(par->eta);

    #pragma omp for schedule(dynamic)
    for(int ib = 0; ib < par->Nb_pd; ib++){

        if(par->max_threads > 1){
            pos_left[1] = 0; // with multiple threads there could be downward jumps in b
        }
        int prev_pos_left_2;
        if(par->ndim == 3){
            prev_pos_left_2 = 0; // k_next is increasing, increasing q loop
        }

    for(int iq = 0; iq < par->Nq_pd; iq++){

        int prev_pos_left_0      = 0;  // m_retire increasing is increasing, increasing a loop
        int prev_pos_left_retire = 0;  // m_next is increasing, increasing a loop

    for(int ia = 0; ia < par->Na_pd; ia++){

        double w  = 0.0;
        double wa = 0.0;
        double wb = 0.0;
        double wq = 0.0;

        pos_left[0]        = prev_pos_left_0;       // from latest i_eta == 0
        pos_left_retire[0] = prev_pos_left_retire;  // from latest i_eta == 0
        if(par->ndim == 3){
            pos_left[2] = prev_pos_left_2;          // from latest i_eta == 0
        }

        // b. loop over shocks
        for(int i_eta = 0; i_eta < par->Neta; i_eta++){

            // i. final next-period states
           double m_next, n_next, k_next, m_retire;
           if(par->ndim == 3){
              m_next    = par->Ra*par->grid_a_pd[ia];
              n_next    = par->Rb*par->grid_b_pd[ib];
              k_next    = par->eta[i_eta]*par->grid_q_pd[iq];
              m_retire  = m_next + n_next + par->rk_retire*k_next;
           } else {
              m_next    = par->Ra*par->grid_a_pd[ia] + par->eta[i_eta];
              n_next    = par->Rb*par->grid_b_pd[ib];
              m_retire  = m_next + n_next;
           }

            // ii. interpolations
            for(int i = 0; i < dimf; i++){
                yi[i] = &(v_vprime_retire[i][i_eta]);
            }
            mesh::interp_multiple_guess_split(v_vprime_retire_interp_func, &m_retire, pos_left_retire, yi);
            next_states[0] = m_next;
            next_states[1] = n_next;
            if(par->ndim == 3){
                next_states[2] = m_next;
            }
            for(int i = 0; i < dimf; i++){
                yi[i] = &(v_vprime[i][i_eta]);
            }
            mesh::interp_multiple_guess_split(v_vprime_interp_func, next_states, pos_left, yi);

                // update position vectors
                if(i_eta == 0){
                    prev_pos_left_0      = pos_left[0];
                    prev_pos_left_retire = pos_left_retire[0];
                    if(par->ndim == 3){
                        prev_pos_left_2  = pos_left[2];
                    }
                }

        } // eta

        // c. vectorized LogSum
        Vec4d v_work_vec, v_retire_vec;

            // load
            v_work_vec.load(v_vprime[0]);
            v_retire_vec.load(v_vprime_retire[0]);

            // trnasform
            v_work_vec   = -1.0/v_work_vec;
            v_retire_vec = -1.0/v_retire_vec;

            // maximum
            Vec4db I = v_work_vec > v_retire_vec;
            Vec4d mxm_vec = select(I,v_work_vec,v_retire_vec);

            Vec4d prob_work_vec, prob_retire_vec, LogSum_vec;
            if(fabs(par->sigma) > 1e-6){

                // contributions
                Vec4d cont_work_vec   = exp((v_work_vec - mxm_vec)*inv_sigma);
                Vec4d cont_retire_vec = exp((v_retire_vec - mxm_vec)*inv_sigma);

                // Logsum
                LogSum_vec  = mxm_vec + par->sigma * log(cont_work_vec + cont_retire_vec);

                // probbilities
                prob_work_vec    = exp((v_work_vec-LogSum_vec)*inv_sigma);
                prob_retire_vec  = 1.0 - prob_work_vec;

            } else {

                LogSum_vec = mxm_vec;

            }

        // d. weighted sum
        w += horizontal_add(w_eta_vec*LogSum_vec);

        if(par->do_derivatives == 1){

            // wa
            Vec4d wa_work_vec, wa_retire_vec;

                // load
                wa_work_vec.load(v_vprime[1]);
                wa_retire_vec.load(v_vprime_retire[1]);

                // transform
                wa_work_vec   = -1.0/wa_work_vec;
                wa_retire_vec = -1.0/wa_retire_vec;

                // sum
                Vec4d wa_vec;
                if(fabs(par->sigma) > 1e-6){
                    wa_vec = prob_work_vec*wa_work_vec + prob_retire_vec*wa_retire_vec;
                } else {
                    wa_vec = select(I,wa_work_vec,wa_retire_vec);
                }
                wa += horizontal_add(w_eta_vec*par->Ra*wa_vec);

            // wb
            Vec4d wb_work_vec, wb_retire_vec;

                // load
                wb_work_vec.load(v_vprime[2]);
                wb_retire_vec.load(v_vprime_retire[2]);

                // work
                wb_work_vec   = -1.0/wb_work_vec;
                wb_retire_vec = -1.0/wb_retire_vec;

                // sum
                Vec4d wb_vec;
                if(fabs(par->sigma) > 1e-6){
                    wb_vec = prob_work_vec*wb_work_vec + prob_retire_vec*wb_retire_vec;
                } else {
                    wb_vec = select(I,wb_work_vec,wb_retire_vec);
                }
                wb += horizontal_add(w_eta_vec*par->Rb*wb_vec);

            // wq
            if(par->ndim == 3){

                Vec4d wq_work_vec, wq_retire_vec;

                    // load
                    wq_work_vec.load(v_vprime[3]);
                    wq_retire_vec.load(v_vprime_retire[3]);

                    // transform
                    wq_work_vec = -1.0/wq_work_vec;
                    wq_retire_vec = -1.0/wq_retire_vec;

                    // sum
                    Vec4d wq_vec;
                    if(fabs(par->sigma) > 1e-6){
                        wq_vec = prob_work_vec*wq_work_vec + prob_retire_vec*wq_retire_vec;
                    } else {
                        wq_vec = select(I,wq_work_vec,wq_retire_vec);
                    }
                    wq += horizontal_add(w_eta_vec * eta_vec * wq_vec);

            }

        }

        // c. save to output
        int igrid = index_func(ia,ib,iq,par->Na_pd,par->Nb_pd,par->Nq_pd);
        w_grid[igrid]  = w;
        wa_grid[igrid] = wa;
        wb_grid[igrid] = wb;
        if(par->ndim == 3){
            wq_grid[igrid] = wq;
        }

    } // b
    } // q
    } // a


        // clean up
        delete[] grids_v_vprime;
        for(int i = 0; i < dimf; i++){
            delete[] v_vprime[i];
        }
        delete[] v_vprime;
        delete[] next_states;
        delete[] grid_vectors;
        delete[] grid_lengths;
        mesh::destroy(v_vprime_interp_func);

        delete[] grid_v_vprime_retire;
        delete[] v_vprime_retire;
        delete[] grid_vector_retire;
        delete[] grid_length_retire;
        mesh::destroy(v_vprime_retire_interp_func);

        delete[] pos_left;
        delete[] pos_left_retire;

    } // parallel

    // clean up
    delete par;
    delete[] dims;

} // gateway
