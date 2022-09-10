/*****************************************************************
 * This programme is used to solve a coupled system describing the 
 * diffusion of chemicals undergoing a reversible reaction 
 * in a rotating chamber with source terms and 'extraction terms'
 * Prior to compilation execute following lines on nenneke:
 * module purge
 * module load intel impi imkl
 * Then:
 * Compile:  gcc -o bandu band_utility.c -lm -lmkl -liomp5
 * Run: ./bandu
 * ****************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <mkl_lapacke.h>

/* Define a new type band_mat */
typedef struct band_mat band_mat;

/* Define structure that holds band matrix information */
struct band_mat{
    long ncol;        /* Number of columns in band matrix */
    long nbrows;      /* Number of rows (bands in original matrix) */
    long nbands_up;   /* Number of bands above diagonal */
    long nbands_low;  /* Number of bands below diagonal */
    double *array;    /* Storage for the matrix in banded format */
    /* Internal temporary storage for solving inverse problem */
    long nbrows_inv;  /* Number of rows of inverse matrix */
    double *array_inv;/* Store the inverse if this is generated */
    int *ipiv;        /* Additional inverse information */
};

int read_coeffs(); 
void read_input(long *N_x, long *N_y, double *L_x, double *L_y, double *t_f, double *t_D, int *s_x0, int *s_x1, int *s_y0, int *s_y1); 
void write_output(double time, double x, double y, double u); 
int init_band_mat(band_mat *bmat, long nbands_lower, long nbands_upper, long n_columns);
void finalise_band_mat(band_mat *bmat);
double *getp(band_mat *bmat, long row, long column);
double getv(band_mat *bmat, long row, long column);
double setv(band_mat *bmat, long row, long column, double val);
int solve_Ax_eq_b(band_mat *bmat, double *x, double *b);
int printmat(band_mat *bmat);
long indx(long k, long m, long K);
int is_valid(long k, long m, long K, long M);
int is_in_bmat(long i, long j, long bands);


int main() {
   /*******************************
    Declare parameters 
    *******************************/
    long N_x; //Number of z grid points
    long N_y; //Number of y grid points
    double L_x; //Length of x domain
    double L_y; //Length of y domain
    double t_f; //final time for time evolution mode
    double t_D; //Diagnostic time step
    int s_x0; //switch for x=0 boundary condition
    int s_x1; //switch for x=L_x boundary condition
    int s_y0; //switch for y=0 boundary condition
    int s_y1; //switch for y=L_y boundary condition

    //Read in values form file
    read_input(&N_x, &N_y, &L_x, &L_y, &t_f, &t_D, &s_x0, &s_x1, &s_y0, &s_y1);

    if((s_x0!=0 && s_x0!=1) || (s_x1!=0 && s_x1!=1) || (s_y0!=0 && s_y0!=1) || (s_y1!=0 && s_y1!=1)){
        printf("Invalid boundary condition switch");
        exit(1);
    }

    double dx = L_x/(N_x-1);
    double dy = L_y/(N_y-1);
    double dt = 0.25/(1/(dx*dx)+1/(dy*dy)); //for stability
    long nsteps = t_D/dt;
    //nsteps+=1;
    dt = t_D/nsteps;


    double alpha;
    double beta;
    double gamma;

    //Allocate memory to sigma and source arrays and check for failed allocation
    double *uc = (double*)malloc(sizeof(double)*N_x*N_y);
    double *ptr1 = uc; //reference pointer
    double *un = (double*)malloc(sizeof(double)*N_x*N_y);
    double *ptr2 = un; //refernce pointer
    

    if(uc==NULL || un==NULL){
        printf("Memory allocation failed");
        return 1;
    }

    //Read from coefficients file
    read_coeffs(uc, N_x, N_y);

    alpha = dt/(2*dx*dx);
    beta = dt/(2*dy*dy);
    gamma = 1+2*alpha+2*beta;

    band_mat bmat;
    long ncols = N_x*N_y;
    
    //Allocate memory to vectors x and b To solve equation Ax=b
    double *x = malloc(sizeof(double)*ncols);
    double *b = malloc(sizeof(double)*ncols);

    if(x==NULL||b==NULL){
        printf("Memory allocation failed");
        return 1;
    }
    
    long k;
    long m;
    long unknown_idx;

    long nbands_up = 2*N_x-1;  
    long nbands_low  = nbands_up;
    init_band_mat(&bmat, nbands_low, nbands_up, ncols);

    /* Default matrix with no boundary conditions */

    for(k=0;k<N_x;k++){
        for(m=0;m<N_y;m++){
            unknown_idx = indx(k, m, N_x);
            if(is_valid(k-1, m, N_x, N_y)){
                setv(&bmat, unknown_idx, indx(k-1, m, N_x), -alpha);
            }
            if(is_valid(k+1, m, N_x, N_y)){
                setv(&bmat, unknown_idx, indx(k+1, m, N_x), -alpha);
            }
            if(is_valid(k, m-1, N_x, N_y)){
                setv(&bmat, unknown_idx, indx(k, m-1, N_x), -beta);
            }
            if(is_valid(k, m+1, N_x, N_y)){
                setv(&bmat, unknown_idx, indx(k, m+1, N_x), -beta);
            }
            setv(&bmat, unknown_idx, indx(k, m, N_x), gamma);

        }
        
    }

    /* BOUNDARY CONDITIONS*/

    // x=0 Dirichlet BC
    long p;
    long set_idx;

    if (s_x0==0){
        for(k=0;k<N_x;k++){
            for(m=0;m<N_y;m++){
                unknown_idx = indx(k, m, N_x);
                for(p=0;p<N_y;p++){
                    set_idx = indx(0, p, N_x);
                    //printf("%ld ", unknown_idx-set_idx+nbands_up);
                    if(is_in_bmat(unknown_idx, set_idx, nbands_up)){
                        setv(&bmat, unknown_idx, set_idx, 0);
                    } 
                }
            }
        }
    }
    // x=0 Neumann BC
    else{
        for(m=0;m<N_y;m++){
            unknown_idx = indx(0, m, N_x);
            setv(&bmat, unknown_idx, indx(1, m, N_x), -2*alpha);
        }
        
    }
    //x=L_x Dirichlet BC
    if (s_x1==0){
        for(k=0;k<N_x;k++){
            for(m=0;m<N_y;m++){
                unknown_idx = indx(k, m, N_x);
                for(p=0;p<N_y;p++){
                    set_idx = indx(N_x-1, p, N_x);
                    //printf("%ld ", unknown_idx-set_idx+nbands_up);
                    if(is_in_bmat(unknown_idx, set_idx, nbands_up)){
                        setv(&bmat, unknown_idx, set_idx, 0);
                    } 
                }
                
            }
        }
    }
    //x=L_x Neumann BC
    else{
        for(m=0;m<N_y;m++){
            unknown_idx = indx(N_x-1, m, N_x);
            setv(&bmat, unknown_idx, indx(N_x-2, m, N_x), -2*alpha);
        }
        
    }

    // y=0 Dirichlet BC
    if(s_y0==0){
        for(k=0;k<N_x;k++){
            for(m=0;m<N_y;m++){
                unknown_idx = indx(k, m, N_x);
                for(p=0;p<N_x;p++){
                    set_idx = indx(p, 0, N_x);
                    //printf("%ld ", unknown_idx-set_idx+nbands_up);
                    if(is_in_bmat(unknown_idx, set_idx, nbands_up)){
                        setv(&bmat, unknown_idx, set_idx, 0);
                    } 
                }
                
            }
        }
    }

    // y=0 Neumann BC
    else{
        for(k=0;k<N_x;k++){
            unknown_idx = indx(k, 0, N_x);
            setv(&bmat, unknown_idx, indx(k, 1, N_x), -2*beta);
        }
    }

    // y=L_y Dirichlet BC
    if(s_y1==0){
        for(k=0;k<N_x;k++){
            for(m=0;m<N_y;m++){
                unknown_idx = indx(k, m, N_x);
                for(p=0;p<N_x;p++){
                    set_idx = indx(p, N_y-1, N_x);
                    //printf("%ld ", unknown_idx-set_idx+nbands_up);
                    if(is_in_bmat(unknown_idx, set_idx, nbands_up)){
                        setv(&bmat, unknown_idx, set_idx, 0);
                    } 
                }
            }
        }
    }

    // y=L_y Neumann BC
    else{
        for(k=0;k<N_x;k++){
            unknown_idx = indx(k, N_y-1, N_x);
            setv(&bmat, k, indx(k, N_y-2, N_x), -2*beta);
        }
    }

    printmat(&bmat);

    double ctime=0;
    long count = 0;
 
        
    while(ctime<t_f + dt/2){       

        for(k=0;k<N_x;k++){
            for(m=0;m<N_y;m++){
                unknown_idx = indx(k, m, N_x);
                b[unknown_idx] = uc[unknown_idx]*(1-2*alpha-2*beta);
                if(is_valid(k-1, m, N_x, N_y)){
                    b[unknown_idx]+=alpha*uc[indx(k-1, m, N_x)];
                }
                if(is_valid(k+1, m, N_x, N_y)){
                    b[unknown_idx]+=alpha*uc[indx(k+1, m, N_x)];
                }
                if(is_valid(k, m-1, N_x, N_y)){
                    b[unknown_idx]+=beta*uc[indx(k, m-1, N_x)];
                }
                if(is_valid(k, m+1, N_x, N_y)){
                    b[unknown_idx]+=beta*uc[indx(k, m+1, N_x)];
                }
                b[unknown_idx]+=(uc[unknown_idx])*(uc[unknown_idx])/(1+(uc[unknown_idx])*(uc[unknown_idx]));

            }
        }
    

        solve_Ax_eq_b(&bmat, x, b);

        //Hardcode Dirichlet BC and update u array

        for(k=0;k<N_x;k++){
            for(m=0;m<N_y;m++){
                unknown_idx = indx(k,m,N_x);
                un[unknown_idx] = x[unknown_idx];
                if(s_x0==0){
                    if(k==0){
                        un[unknown_idx]=0;
                    }
                }
                if(s_x1==0){
                    if(k==N_x-1){
                        un[unknown_idx]=0;
                    }
                }
                if(s_y0==0){
                    if(m==0){
                        un[unknown_idx]=0;
                    }
                }
                if(s_y1==0){
                    if(m==N_y-1){
                        un[unknown_idx]=0;
                    }
                }
                
            }
        }

        //Output results when a multiple of diagnostic time step
        if(count%nsteps==0){
            for(m=0;m<N_y;m++){
                for(k=0;k<N_x;k++){
                    unknown_idx = indx(k, m, N_x);
                    write_output(ctime, k*dx, m*dy, uc[unknown_idx]);
                    //printf("%lf %lf %lf %lf \n", ctime, k*dx, m*dy, uc[unknown_idx]);
                    uc[unknown_idx] = un[unknown_idx];
                }
            }
        }

        ctime+=dt;
        count+=1;
        
    }

    uc = ptr1;
    un = ptr2;

    
    finalise_band_mat(&bmat);
    free(x);
    free(b);
    free(uc);
    free(un);

    return(0);
}

/* k is x coordinate
    m is y coordinate
    K is number of x grid points
*/

long indx(long k, long m, long K){
    return m*K + k;
}

/*Check that a grid point has valid coordinates */
int is_valid(long k, long m, long K, long M) {
  return (k>=0)&&(k<K)&&(m>=0)&&(m<M);
}

//Checking that the value is in inside re-indexed matrix so no unnecessary appends
int is_in_bmat(long i, long j, long bands){
    return (i-j+bands>=0)&&(i-j<=bands);
}

int read_coeffs(double *uc, long N_x, long N_y) {

    FILE *coeff_file;
    if(!(coeff_file=fopen("coefficients.txt","r"))){
        printf("Error opening file\n");
        exit(1);
    }

    long j;
    for(j=0;j<N_x*N_y;j++) {
        fscanf(coeff_file,"%lf", &uc[j]);
    }
    fclose(coeff_file);

    return 0;
}

void read_input(long *N_x, long *N_y, double *L_x, double *L_y, double *t_f, double *t_D, int *s_x0, int *s_x1, int *s_y0, int *s_y1) {
    FILE *infile;
    if(!(infile=fopen("input.txt","r"))) {
        printf("Error opening file\n");
        exit(1);
    }
    if(10!=fscanf(infile,"%ld %ld %lf %lf %lf %lf %d %d %d %d", N_x, N_y, L_x, L_y, t_f, t_D, s_x0, s_x1, s_y0, s_y1)) {
        //Entries are on different lines so may need to edit
        printf("Error reading parameters from file\n");
        exit(1);
    }
    fclose(infile);
}

void write_output(double time, double x, double y, double u){
    FILE *outfile;
    if(!(outfile=fopen("output.txt", "a"))){
        printf("Error opening output file\n");
        exit(1);
    }
    fprintf(outfile, "%lf %lf %lf %lf\n", time, x, y, u);
    fclose(outfile);
}

/* Initialise a band matrix of a certain size, allocate memory,
   and set the parameters.  */ 
int init_band_mat(band_mat *bmat, long nbands_lower, long nbands_upper, long n_columns) {
    bmat->nbrows = nbands_lower + nbands_upper + 1;
    bmat->ncol   = n_columns;
    bmat->nbands_up = nbands_upper;
    bmat->nbands_low= nbands_lower;
    bmat->array      = (double *) malloc(sizeof(double)*bmat->nbrows*bmat->ncol);
    bmat->nbrows_inv = bmat->nbands_up*2 + bmat->nbands_low + 1;
    bmat->array_inv  = (double *) malloc(sizeof(double)*(bmat->nbrows+bmat->nbands_low)*bmat->ncol);
    bmat->ipiv       = (int *) malloc(sizeof(int)*bmat->ncol);
    if (bmat->array==NULL||bmat->array_inv==NULL) {
        return 0;
    }  
  /* Initialise array to zero */
    long i;
    for (i=0;i<bmat->nbrows*bmat->ncol;i++) {
        bmat->array[i] = 0.0;
    }
    return 1;
}

/* Finalise function: should free memory as required */
void finalise_band_mat(band_mat *bmat) {
    free(bmat->array);
    free(bmat->array_inv);
    free(bmat->ipiv);
}

/* Get a pointer to a location in the band matrix, using
   the row and column indexes of the full matrix.           */
double *getp(band_mat *bmat, long row, long column) {
    int bandno = bmat->nbands_up + row - column;
    if(row<0 || column<0 || row>=bmat->ncol || column>=bmat->ncol ) {
        printf("Indexes out of bounds in getp: %ld %ld %ld \n",row,column,bmat->ncol);
        exit(1);
    }
    return &bmat->array[bmat->nbrows*column + bandno];
}

/* Retrun the value of a location in the band matrix, using
   the row and column indexes of the full matrix.           */
double getv(band_mat *bmat, long row, long column) {
    return *getp(bmat,row,column);
}

/* Set an element of a band matrix to a desired value based on the pointer
   to a location in the band matrix, using the row and column indexes
   of the full matrix.           */
double setv(band_mat *bmat, long row, long column, double val) {
    *getp(bmat,row,column) = val;
    return val;
}

/* Solve the equation Ax = b for a matrix a stored in band format
   and x and b real arrays                                          */
int solve_Ax_eq_b(band_mat *bmat, double *x, double *b) {
  /* Copy bmat array into the temporary store */
    int i,bandno;
    for(i=0;i<bmat->ncol;i++) { 
        for (bandno=0;bandno<bmat->nbrows;bandno++) {
            bmat->array_inv[bmat->nbrows_inv*i+(bandno+bmat->nbands_low)] = bmat->array[bmat->nbrows*i+bandno];
        }
        x[i] = b[i];
    }
    long nrhs = 1;
    long ldab = bmat->nbands_low*2 + bmat->nbands_up + 1;
    int info = LAPACKE_dgbsv( LAPACK_COL_MAJOR, bmat->ncol, bmat->nbands_low, bmat->nbands_up, nrhs, bmat->array_inv, ldab, bmat->ipiv, x, bmat->ncol);
    return info;
}

int printmat(band_mat *bmat) {
    long i,j;
    for(i=0; i<bmat->ncol;i++) {
        for(j=0; j<bmat->nbrows; j++) {
            printf("%ld %ld %g \n",i,j,bmat->array[bmat->nbrows*i + j]);
        }
    }
    return 0;
}



