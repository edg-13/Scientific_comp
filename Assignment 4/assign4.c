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

int read_coeffs(); //read in from coefficients.txt
void read_input(double *L, long *N, double *D, double *v, double *k_plus, double *k_minus); //read in from input.txt
int init_band_mat(band_mat *bmat, long nbands_lower, long nbands_upper, long n_columns);
void finalise_band_mat(band_mat *bmat);
double *getp(band_mat *bmat, long row, long column);
double getv(band_mat *bmat, long row, long column);
double setv(band_mat *bmat, long row, long column, double val);
int solve_Ax_eq_b(band_mat *bmat, double *x, double *b);
int printmat(band_mat *bmat);


int main() {
   /*******************************
    Declare parameters 
    *******************************/
    double L; //Right boundary of x points
    long N; //Number of grid points
    double D; //Diffusion coefficient
    double v; //Advection velocity
    double k_plus; //Forward rate constant of reaction
    double k_minus; //Backward rate constant of reaction

    double alpha1; //a_k-1 terms in equation for A and B for v>0
    double alpha2; //a_k-1 terms in equation for A and B for v<0
    double beta_a1; //a_k terms in equation for A with v>0
    double beta_b1; //a_k terms in equation for B with v>0 without sigma
    double beta_a2; //a_k terms in equation for A with v<0
    double beta_b2; //a_k terms in equation for B with v<0 without sigma
    double gamma; //a_k+1 terms in equations A and B

    //Read in values form file
    read_input(&L, &N, &D, &v, &k_plus, &k_minus);

    //Allocate memory to sigma and source arrays and check for failed allocation
    double *S = (double*)malloc(sizeof(double)*N);
    double *sigma = (double*)malloc(sizeof(double)*N);

    if(S==NULL || sigma==NULL){
        printf("Memory allocation failed");
        return 1;
    }

    //Read from sigma and S values from coefficients file
    read_coeffs(S, sigma, N);

    //x domain is x_0 to x_(N-1)
    double dx = L/N;
    double inv_dx = 1.0/dx;

    band_mat bmat;
    long ncols = N*2;
    
    long nbands_low = 4;  
    long nbands_up  = 4;
    init_band_mat(&bmat, nbands_low, nbands_up, ncols);

    //Allocate memory to vectors x and b To solve equation Ax=b
    double *x = malloc(sizeof(double)*ncols);
    double *b = malloc(sizeof(double)*ncols);

    if(x==NULL||b==NULL){
        printf("Memory allocation failed");
        return 1;
    }

    //Advective term switches sign in upwinding scheme depending on sign of v
    alpha1 = D*inv_dx*inv_dx + v*inv_dx; 
    alpha2 = D*inv_dx*inv_dx - v*inv_dx; 
    beta_a1 = -2*D*inv_dx*inv_dx - v*inv_dx - k_plus; 
    beta_b1 = -2*D*inv_dx*inv_dx - v*inv_dx - k_minus; 
    beta_a2 = -2*D*inv_dx*inv_dx + v*inv_dx - k_plus; 
    beta_b2 = -2*D*inv_dx*inv_dx + v*inv_dx - k_minus; 
    gamma = D*inv_dx*inv_dx; 

    long i;
    /* Loop over the equation number and set the matrix
        values equal to the coefficients of the grid values. 
        Note that boundaries are treated with special cases */
    
    //Matrix is set up as FOLLOWS: indices are arranged as a_0, b_0, a_(N-1), b_(N-1), a_1, b_1, ...


    if(v>0){
        for(i=0; i<ncols; i++) {
            if(i%2==0){
                setv(&bmat,i,i,beta_a1);
            }
            else{
                if((i-1)%4==0){
                    setv(&bmat,i,i,beta_b1 - sigma[(i-1)/4]);
                }
                else{
                    setv(&bmat,i,i,beta_b1 - sigma[N-1-((i-3)/4)]);
                }
            }
            if(i>0){
                if(i%2==1){
                    setv(&bmat,i,i-1,k_plus);
                }
                else{
                    setv(&bmat,i,i-1,0);
                }
            }
            if(i>1){
                setv(&bmat,i,i-2,0);
            }
            if(i>2){
                setv(&bmat,i,i-3,0);
            }
            if(i>3){
                if(i%4==0){
                    setv(&bmat,i,i-4,alpha1);
                    setv(&bmat,i+1,i-3,alpha1);
                }
                else if((i+2)%4==0){
                    setv(&bmat,i,i-4,gamma);
                    setv(&bmat,i+1,i-3,gamma);
                }
            }
            if(i<ncols-1){
                if(i%2==0){
                    setv(&bmat,i,i+1,k_minus);
                }
                else{
                    setv(&bmat,i,i+1,0);
                }
            }
            if(i<ncols-2){
                setv(&bmat,i,i+2,0);
            }
            if(i<ncols-3){
                setv(&bmat,i,i+3,0);
            }
            if(i<ncols-4){
                if(i%4==0){
                    setv(&bmat,i,i+4,gamma);
                    setv(&bmat,i+1,i+5,gamma);
                }
                else if((i+2)%4==0){
                    setv(&bmat,i,i+4,alpha1);
                    setv(&bmat,i+1,i+5,alpha1);
                }
            }
            if(i%2==1){
                b[i]=0;
            }
            else{
                if(i%4==0){
                    b[i]=-S[i/4];
                }
                else{
                    b[i]=-S[N-1-((i-2)/4)];
                }
            }
        }
        setv(&bmat,2,0,gamma);
        setv(&bmat,3,1,gamma);

        setv(&bmat,0,2,alpha1);
        setv(&bmat,1,3,alpha1);

        if(ncols%4==0){
            setv(&bmat,ncols-1,ncols-3,alpha1);
            setv(&bmat,ncols-2,ncols-4,alpha1);

            setv(&bmat,ncols-3,ncols-1,gamma);
            setv(&bmat,ncols-4,ncols-2,gamma);
        }
        else{
            setv(&bmat,ncols-1,ncols-3,gamma);
            setv(&bmat,ncols-2,ncols-4,gamma);

            setv(&bmat,ncols-3,ncols-1,alpha1);
            setv(&bmat,ncols-4,ncols-2,alpha1);
        }
    }   

    else{
        for(i=0; i<ncols; i++) {
            if(i%2==0){
                setv(&bmat,i,i,beta_a2);
            }
            else{
                if((i-1)%4==0){
                    setv(&bmat,i,i,beta_b2 - sigma[(i-1)/4]);
                }
                else{
                    setv(&bmat,i,i,beta_b2 - sigma[N-1-((i-3)/4)]);
                }
            }
            if(i>0){
                if(i%2==1){
                    setv(&bmat,i,i-1,k_plus);
                }
                else{
                    setv(&bmat,i,i-1,0);
                }
            }
            if(i>1){
                setv(&bmat,i,i-2,0);
            }
            if(i>2){
                setv(&bmat,i,i-3,0);
            }
            if(i>3){
                if(i%4==0){
                    setv(&bmat,i,i-4,alpha2);
                    setv(&bmat,i+1,i-3,alpha2);
                }
                else if((i+2)%4==0){
                    setv(&bmat,i,i-4,gamma);
                    setv(&bmat,i+1,i-3,gamma);
                }
            }
            if(i<ncols-1){
                if(i%2==0){
                    setv(&bmat,i,i+1,k_minus);
                }
                else{
                    setv(&bmat,i,i+1,0);
                }
            }
            if(i<ncols-2){
                setv(&bmat,i,i+2,0);
            }
            if(i<ncols-3){
                setv(&bmat,i,i+3,0);
            }
            if(i<ncols-4){
                if(i%4==0){
                    setv(&bmat,i,i+4,gamma);
                    setv(&bmat,i+1,i+5,gamma);
                }
                else if((i+2)%4==0){
                    setv(&bmat,i,i+4,alpha2);
                    setv(&bmat,i+1,i+5,alpha2);
                }
            }
            if(i%2==1){
                b[i]=0;
            }
            else{
                if(i%4==0){
                    b[i]=-S[i/4];
                }
                else{
                    b[i]=-S[N-1-((i-2)/4)];
                }
            }
        }
        setv(&bmat,2,0,gamma);
        setv(&bmat,3,1,gamma);

        setv(&bmat,0,2,alpha2);
        setv(&bmat,1,3,alpha2);

        if(ncols%4==0){
            setv(&bmat,ncols-1,ncols-3,alpha2);
            setv(&bmat,ncols-2,ncols-4,alpha2);

            setv(&bmat,ncols-3,ncols-1,gamma);
            setv(&bmat,ncols-4,ncols-2,gamma);
        }
        else{
            setv(&bmat,ncols-1,ncols-3,gamma);
            setv(&bmat,ncols-2,ncols-4,gamma);

            setv(&bmat,ncols-3,ncols-1,alpha2);
            setv(&bmat,ncols-4,ncols-2,alpha2);
        }
    }
  
    /*  Print coefficient matrix for debugging: */ 
    //printmat(&bmat);

    //long k;
    double val1;
    double val2;    

    solve_Ax_eq_b(&bmat, x, b);


    //'Unfolding' the x domain to retrieve original ordering of indices
    for(i=0;i<N;i++){
        if(N%2==0){
            if(i<N/2){
                val1 = x[4*i];
                val2 = x[(4*i)+1];
            }
            else{
                val1 = x[4*(N-i)-2];
                val2 = x[4*(N-i)-1];
            }
        }
        else{
            if(i<=N/2){
                val1 = x[4*i];
                val2 = x[(4*i)+1];
            }
            else{
                val1 = x[4*(N-i)-2];
                val2 = x[4*(N-i)-1];
            }

        }
        printf("%lf %lf %lf\n", i*dx, val1, val2);
    }

    /***********************
    Freeing allocated memory
    ***********************/
    finalise_band_mat(&bmat); 
    free(S);
    free(sigma);
    free(x);
    free(b);
    return(0);
}

int read_coeffs(double *S_arr, double *sigma_arr, long N) {

    FILE *coeff_file;
    if(!(coeff_file=fopen("coefficients.txt","r"))){
        printf("Error opening file\n");
        exit(1);
    }

    long j;
    for(j=0;j<N;j++) {
        fscanf(coeff_file,"%lf      %lf", &S_arr[j], &sigma_arr[j]);
    }
    fclose(coeff_file);

    return 0;
}

void read_input(double *L, long *N, double *D, double *v, double *k_plus, double *k_minus) {
    FILE *infile;
    if(!(infile=fopen("input.txt","r"))) {
        printf("Error opening file\n");
        exit(1);
    }
    if(6!=fscanf(infile,"%lf %ld %lf %lf %lf %lf", L, N, D, v, k_plus, k_minus)) {
        printf("Error reading parameters from file\n");
        exit(1);
    }
    fclose(infile);
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



