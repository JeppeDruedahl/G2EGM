#include<cmath>

namespace mesh {

typedef struct settings_struct
{
    // input
    int dim, dimf, *Nx;
    double **x, *y, **y_multiple;

    // pre-calculations
    int ncube;
    int *pos_left, **add, *facs;
    double *reldiff, **cs;

} settings_struct;

void create(mesh::settings_struct *settings, int dim, int *Nx, double **x, double *y);
void create_multiple(mesh::settings_struct *settings, int dim, int *Nx, double **x, int dimf, double **y);
void destroy(mesh::settings_struct *settings);
int binary_search(int imin, int Nx, double *x, double xi);
double interp(mesh::settings_struct *settings, double *xi);
void interp_multiple(mesh::settings_struct *settings, double *xi, double *yi);
void interp_multiple_guess(mesh::settings_struct *settings, double *xi, int *pos_left_guess, double *yi);
void interp_multiple_guess_split(mesh::settings_struct *settings, double *xi, int *pos_left_guess, double **yi);

}

void mesh::create(mesh::settings_struct *settings, int dim, int *Nx, double **x, double *y)
{

    int i, j, k, *add;

    // a. basics
    settings->dim         = dim;
    settings->ncube       = (int)pow(2.0,(double)(dim-1));
    settings->Nx          = Nx;
    settings->x           = x;
    settings->y           = y;

    // b. position factors
    settings->facs = new int[dim];
    for(i = 0; i < dim; i++){
        settings->facs[i] = 1;
        for(j = 0; j < i; j++){
            settings->facs[i] *= Nx[j];
        }
    }

    // c. containers
    settings->pos_left  = new int[dim];
    settings->reldiff   = new double[dim];
    settings->cs        = new double*[dim];
    for(j = 0; j < dim; j++){
        settings->cs[j] = new double[(int)pow(2.0,(double)(dim-1-j))];
    }

    // d. add (ordering of hypercube corners is important)
    settings->add = new int*[settings->ncube];
    for(i = 0; i < settings->ncube; i++){

        settings->add[i] = new int[dim];
        add          = settings->add[i];

        if(i == 0){
            for(j = 0; j < dim; j++){ add[j] = 0;}
        } else {
            for(j = 0; j < dim; j++){ add[j] = settings->add[i-1][j]; }
            for(j = 1; j < dim; j++){
                if(add[j] == 0){
                    add[j] = 1;
                    for(k = 1; k < j; k++){ add[k] = 0; }
                    break;
                }
            }
        } // if i == 0

    } // i

} // create

void mesh::create_multiple(mesh::settings_struct *settings, int dim, int *Nx, double **x, int dimf, double **y_multiple)
{

    int i, j, k, *add;

    // a. basics
    settings->dim         = dim;
    settings->dimf        = dimf;
    settings->ncube       = (int)pow(2.0,(double)(dim-1));
    settings->Nx          = Nx;
    settings->x           = x;
    settings->y_multiple  = y_multiple;

    // b. position factors
    settings->facs = new int[dim];
    for(i = 0; i < dim; i++){
        settings->facs[i] = 1;
        for(j = 0; j < i; j++){
            settings->facs[i] *= Nx[j];
        }
    }

    // c. containers
    settings->pos_left  = new int[dim];
    settings->reldiff   = new double[dim];
    settings->cs        = new double*[dim];
    for(j = 0; j < dim; j++){
        settings->cs[j] = new double[dimf*(int)pow(2.0,(double)(dim-1-j))];
    }

    // d. add (ordering of hypercube corners is important)
    settings->add = new int*[settings->ncube];
    for(i = 0; i < settings->ncube; i++){

        settings->add[i] = new int[dim];
        add          = settings->add[i];

        if(i == 0){
            for(j = 0; j < dim; j++){ add[j] = 0;}
        } else {
            for(j = 0; j < dim; j++){ add[j] = settings->add[i-1][j]; }
            for(j = 1; j < dim; j++){
                if(add[j] == 0){
                    add[j] = 1;
                    for(k = 1; k < j; k++){ add[k] = 0; }
                    break;
                }
            }
        } // if i == 0

    } // i

} // create

void mesh::destroy(mesh::settings_struct *settings)
{

    int i, d, dim = settings->dim;

    // position factors
    delete[] settings->facs;

    // containers
    delete[] settings->pos_left;
    delete[] settings->reldiff;
    for(d = 0; d < dim; d++){
        delete[] settings->cs[d];;
    }
    delete[] settings->cs;

    // add
    for(i = 0; i < settings->ncube; i++){
        delete[] settings->add[i];
    }
    delete[] settings->add;

    // itself
    delete settings;

} // delete_settings

int mesh::binary_search(int imin, int Nx, double *x, double xi)
{

    int imid, half;

    // checks
    if(xi <= x[0]){
        return 0;
    } else if(xi >= x[Nx-2]) {
        return Nx-2;
    }

    // binary search
    while((half = Nx/2)){
        imid = imin + half;
        imin = (x[imid] <= xi) ? imid:imin;
        Nx  -= half;
    }

    return imin;


} // binary search

double mesh::interp(mesh::settings_struct *settings, double *xi)
{

    int d, i, j;
    int dim, *Nx, *pos_left, index_left, n;
    double **x, *y, *xvec, *reldiff, *c, *cprev;

        ////////////////////////
        // load from settings //
        ////////////////////////

        dim = settings->dim;
        Nx  = settings->Nx;
        x   = settings->x;
        y   = settings->y;

        pos_left  = settings->pos_left;
        reldiff   = settings->reldiff;

    ////////////////////////////////////////////
    // 1. relative distance in each dimension //
    ////////////////////////////////////////////

    for(d = 0; d < dim; d++){

        // a. x-vector in dimension d
        xvec = x[d];

        // b. position to the left of xi in dimension d
        pos_left[d] = mesh::binary_search(0, Nx[d], xvec, xi[d]);

        // c. relative position if xi between neighboring points
        reldiff[d] = (xi[d] - xvec[pos_left[d]]) / (xvec[pos_left[d]+1] - xvec[pos_left[d]]);

    }

    //////////////////////////////////////////////////
    // 2. find function values at hypercube corners //
    //////////////////////////////////////////////////

    n   = settings->ncube;
    c   = settings->cs[0];

    d = 0;
    for(i = 0; i < n; i++){

        index_left = 0;
        for(j = 0; j < dim; j++){ // vectorized by the compiler
            index_left += settings->facs[j]*(pos_left[j] + settings->add[i][j]);
        }


        c[i] = y[index_left]*(1.0-reldiff[d])+y[index_left+1]*reldiff[d];

    }

    ///////////////////////////////////////
    // 3. sequentially remove dimensions //
    ///////////////////////////////////////

    for(d = 1; d < dim; d++){

        // a. number of corners
        n /= 2;

        // b. calculate c in current dimension
        cprev = c;
        c     = settings->cs[d];
        for(i = 0; i < n; i++){ // vectorized by the compiler
            c[i] = cprev[2*i]*(1.0-reldiff[d])+cprev[2*i+1]*reldiff[d];
        }

    }

    // return zero dimensional result
    return c[0];

}

void mesh::interp_multiple(mesh::settings_struct *settings, double *xi, double *yi)
{

    int d, i, j, k;
    int dim, dimf, *Nx, *pos_left, index_left, n, nprev;
    double **x, **y, *xvec, *reldiff, *c, *cprev;

        ////////////////////////
        // load from settings //
        ////////////////////////

        dim  = settings->dim;
        dimf = settings->dimf;
        Nx   = settings->Nx;
        x    = settings->x;
        y    = settings->y_multiple;

        pos_left  = settings->pos_left;
        reldiff   = settings->reldiff;

    ////////////////////////////////////////////
    // 1. relative distance in each dimension //
    ////////////////////////////////////////////



    for(d = 0; d < dim; d++){

        // a. x-vector in dimension d
        xvec = x[d];

        // b. position to the left of xi in dimension d
        pos_left[d] = mesh::binary_search(0, Nx[d], xvec, xi[d]);

        // c. relative position if xi between neighboring points
        reldiff[d] = (xi[d] - xvec[pos_left[d]]) / (xvec[pos_left[d]+1] - xvec[pos_left[d]]);

    }


    //////////////////////////////////////////////////
    // 2. find function values at hypercube corners //
    //////////////////////////////////////////////////

    d   = 0;
    n   = settings->ncube;
    c   = settings->cs[0];

    for(i = 0; i < n; i++){

        index_left = 0;
        for(j = 0; j < dim; j++){ // vectorized by the compiler
            index_left += settings->facs[j]*(pos_left[j] + settings->add[i][j]);
        }


        for(k = 0; k < dimf; k++){
            c[k*n + i] = y[k][index_left]*(1.0-reldiff[d])+y[k][index_left+1]*reldiff[d];
        }

    }

    ///////////////////////////////////////
    // 3. sequentially remove dimensions //
    ///////////////////////////////////////

    for(d = 1; d < dim; d++){

        // a. number of corners
        nprev = n;
        n /= 2;

        // b. calculate c in current dimension
        cprev = c;
        c     = settings->cs[d];
        for(k = 0; k < dimf; k++){
        for(i = 0; i < n; i++){ // vectorized by the compiler
            c[k*n + i] = cprev[k*nprev + 2*i]*(1.0-reldiff[d])+cprev[k*nprev + 2*i+1]*reldiff[d];
        }
        }

    }

    // return zero dimensional result
    for(k = 0; k < dimf; k++){
        yi[k] = c[k];
    }

}

void mesh::interp_multiple_guess(mesh::settings_struct *settings, double *xi, int *pos_left_guess, double *yi)
{

    int d, i, j, k;
    int dim, dimf, *Nx, *pos_left, index_left, n, nprev;
    double **x, **y, *xvec, *reldiff, *c, *cprev;

        ////////////////////////
        // load from settings //
        ////////////////////////

        dim  = settings->dim;
        dimf = settings->dimf;
        Nx   = settings->Nx;
        x    = settings->x;
        y    = settings->y_multiple;

        pos_left  = settings->pos_left;
        reldiff   = settings->reldiff;

    ////////////////////////////////////////////
    // 1. relative distance in each dimension //
    ////////////////////////////////////////////

    for(d = 0; d < dim; d++){

        // a. x-vector in dimension d
        xvec = x[d];

        // b. position to the left of xi in dimension d

            // i. zero guess: binary search
            if(pos_left_guess[d] == 0){
                pos_left[d]       = mesh::binary_search(0, Nx[d], xvec, xi[d]);
            }
            // ii. non-zero guess: while search
            else {
                i = pos_left_guess[d];
                while(xi[d] > xvec[i+1] && i < Nx[d]-2){ i++;}
                pos_left[d] = i;
            }

            // iii. saving the updated guess
            pos_left_guess[d] = pos_left[d];

        // c. relative position if xi between neighboring points
        reldiff[d] = (xi[d] - xvec[pos_left[d]]) / (xvec[pos_left[d]+1] - xvec[pos_left[d]]);

    }


    //////////////////////////////////////////////////
    // 2. find function values at hypercube corners //
    //////////////////////////////////////////////////

    d   = 0;
    n   = settings->ncube;
    c   = settings->cs[0];

    for(i = 0; i < n; i++){

        index_left = 0;
        for(j = 0; j < dim; j++){ // vectorized by the compiler
            index_left += settings->facs[j]*(pos_left[j] + settings->add[i][j]);
        }


        for(k = 0; k < dimf; k++){
            c[k*n + i] = y[k][index_left]*(1.0-reldiff[d])+y[k][index_left+1]*reldiff[d];
        }

    }

    ///////////////////////////////////////
    // 3. sequentially remove dimensions //
    ///////////////////////////////////////

    for(d = 1; d < dim; d++){

        // a. number of corners
        nprev = n;
        n /= 2;

        // b. calculate c in current dimension
        cprev = c;
        c     = settings->cs[d];
        for(k = 0; k < dimf; k++){
        for(i = 0; i < n; i++){ // vectorized by the compiler
            c[k*n + i] = cprev[k*nprev + 2*i]*(1.0-reldiff[d])+cprev[k*nprev + 2*i+1]*reldiff[d];
        }
        }

    }

    // return zero dimensional result
    for(k = 0; k < dimf; k++){
        yi[k] = c[k];
    }

}

void mesh::interp_multiple_guess_split(mesh::settings_struct *settings, double *xi, int *pos_left_guess, double **yi)
{

    int d, i, j, k;
    int dim, dimf, *Nx, *pos_left, index_left, n, nprev;
    double **x, **y, *xvec, *reldiff, *c, *cprev;

        ////////////////////////
        // load from settings //
        ////////////////////////

        dim  = settings->dim;
        dimf = settings->dimf;
        Nx   = settings->Nx;
        x    = settings->x;
        y    = settings->y_multiple;

        pos_left  = settings->pos_left;
        reldiff   = settings->reldiff;

    ////////////////////////////////////////////
    // 1. relative distance in each dimension //
    ////////////////////////////////////////////

    for(d = 0; d < dim; d++){

        // a. x-vector in dimension d
        xvec = x[d];

        // b. position to the left of xi in dimension d

            // i. zero guess: binary search
            if(pos_left_guess[d] == 0){
                pos_left[d]       = mesh::binary_search(0, Nx[d], xvec, xi[d]);
            }
            // ii. non-zero guess: while search
            else {
                i = pos_left_guess[d];
                while(xi[d] > xvec[i+1] && i < Nx[d]-2){ i++;}
                pos_left[d] = i;
            }

            // iii. saving the updated guess
            pos_left_guess[d] = pos_left[d];

        // c. relative position if xi between neighboring points
        reldiff[d] = (xi[d] - xvec[pos_left[d]]) / (xvec[pos_left[d]+1] - xvec[pos_left[d]]);

    }


    //////////////////////////////////////////////////
    // 2. find function values at hypercube corners //
    //////////////////////////////////////////////////

    d   = 0;
    n   = settings->ncube;
    c   = settings->cs[0];

    for(i = 0; i < n; i++){

        index_left = 0;
        for(j = 0; j < dim; j++){ // vectorized by the compiler
            index_left += settings->facs[j]*(pos_left[j] + settings->add[i][j]);
        }


        for(k = 0; k < dimf; k++){
            c[k*n + i] = y[k][index_left]*(1.0-reldiff[d])+y[k][index_left+1]*reldiff[d];
        }

    }

    ///////////////////////////////////////
    // 3. sequentially remove dimensions //
    ///////////////////////////////////////

    for(d = 1; d < dim; d++){

        // a. number of corners
        nprev = n;
        n /= 2;

        // b. calculate c in current dimension
        cprev = c;
        c     = settings->cs[d];
        for(k = 0; k < dimf; k++){
        for(i = 0; i < n; i++){ // vectorized by the compiler
            c[k*n + i] = cprev[k*nprev + 2*i]*(1.0-reldiff[d])+cprev[k*nprev + 2*i+1]*reldiff[d];
        }
        }

    }

    // return zero dimensional result
    for(k = 0; k < dimf; k++){
        yi[k][0] = c[k];
    }

}
