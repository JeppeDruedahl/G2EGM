#define MAX(X,Y) ((X)>(Y)?(X):(Y))
#define MIN(X,Y) ((X)<(Y)?(X):(Y))

/////////////
// STRUCTS //
/////////////

typedef struct
{

    // grids
    int ndim;
    double *grid_m, *grid_n, *grid_k;
    int Nm, Nn, Nk;
    double *grid_a_pd, *grid_b_pd, *grid_q_pd;
    int Na_pd, Nb_pd, Nq_pd;
    int Nm_ret, Neta, N_guess_vfi;

    // parameters
    double beta, rho, alpha, varphi, gamma, chi;
    double rk, rk_retire, delta, lmax;
    double Ra, Rb, yret, sigma;
    double *eta, *w_eta;

    // settings
    int do_derivatives;
    int max_threads, egm_extrap_add;
    double egm_extrap_w;

} par_struct;


////////////////////
// INDEX FUNCTION //
////////////////////

int index_func(int im, int in, int ik, int Nm, int Nn, int Nk)
{
    return ik*Nm*Nn + in*Nm + im;
}

//////////////////////
// UTILITY FUNCTION //
//////////////////////

double u(double c, double l, double rho, double alpha, double varphi, double gamma)
{
    return pow(c , 1.0-rho)/(1.0-rho) - alpha - varphi*pow(l, 1.0+gamma)/(1.0+gamma);
}

////////////
// LOGSUM //
////////////

void logsum(double *LogSum, double *prob1, double *prob2, double *v1, double *v2, double sigma, double inv_sigma, int N)
{

    int i;
    double mxm, cont1, cont2;

    for(i = 0; i < N; i++){

        // a. find maximum
        mxm = MAX(v1[i],v2[i]);

        // b. LogSum and choice probabilities
        if(fabs(sigma) > 0){

            cont1 = exp((v1[i] - mxm)*inv_sigma);
            cont2 = exp((v2[i] - mxm)*inv_sigma);

            LogSum[i] = mxm + sigma*log(cont1 + cont2);

            prob1[i] = exp((v1[i]-LogSum[i])*inv_sigma);
            prob2[i] = exp((v2[i]-LogSum[i])*inv_sigma);

        } else {

            LogSum[i] = mxm;
            if(v1[i] > v2[i]){
                prob1[i]  = 1.0;
                prob2[i]  = 0.0;
            } else {
                prob1[i]  = 0.0;
                prob2[i]  = 1.0;
            }

        }

    } // loop

} // logsum
