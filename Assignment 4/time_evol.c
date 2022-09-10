#include <stdlib.h>
#include <stdio.h>


int read_coeffs();
void read_input(double *L, int *N, double *D, double *v, double *k_plus, double *k_minus);


int main(){

    /*Declare parameters */
    double L; //Right boundary of x points
    int N; //Number of grid points
    double D; //Diffusion coefficient
    double v; //Advection velocity
    double k_plus; //Forward rate constant of reaction
    double k_minus; //Backward rate constant of reaction
   
    read_input(&L, &N, &D, &v, &k_plus, &k_minus);
    
    double dx = L/N;
    double dt = 0.01;
    double final = 60.0;
    
    double *S = (double*)malloc(sizeof(double)*N);
    double *sigma = (double*)malloc(sizeof(double)*N);

    read_coeffs(S, sigma, N);

    double *ac, *an, *bc, *bn, *ac_ptr, *bc_ptr;

    ac = (double *)(malloc(sizeof(double)*N));
    ac_ptr = ac;
    an = (double *)(malloc(sizeof(double)*N));
    bc = (double *)(malloc(sizeof(double)*N));
    bc_ptr = bc;
    bn = (double *)(malloc(sizeof(double)*N));

    if (ac==NULL||an==NULL||bc==NULL||bn==NULL) {
        printf("Memory allocation failed\n");
        return 1;
    }

    int k;
    double x;
    /* Current time */
    double ctime = 0;
    double deriv1;
    double deriv2;
    /* Initialise arrays */
    for(k = 0; k < N; k++) {
        x = k*dx;
        ac[k]  = 3.0;
        bc[k]  = 1.0;
    }
    //printf("#");

    while(ctime < final + dt/2){
        for(k=0; k<N; k++){
            int km = (k+1)%N;
            int kp = (k+N-1)%N;
            if(v>0){
                deriv1 = D*(ac[km]+ac[kp]-2*ac[k])/(dx*dx) - v*(ac[k]-ac[kp])/dx;
                deriv2 = D*(bc[km]+bc[kp]-2*bc[k])/(dx*dx) - v*(bc[k]-bc[kp])/dx;
            }
            else{
                deriv1 = D*(ac[km]+ac[kp]-2*ac[k])/(dx*dx) + v*(ac[k]-ac[kp])/dx;
                deriv2 = D*(bc[km]+bc[kp]-2*bc[k])/(dx*dx) + v*(bc[k]-bc[kp])/dx;
            }
            an[k] = ac[k] + (deriv1 - k_plus*ac[k] + k_minus*bc[k] + S[k])*dt;
            bn[k] = bc[k] + (deriv2 + k_plus*ac[k] - k_minus*bc[k] - sigma[k]*bc[k])*dt;
        }
        ctime+=dt;
        for(k=0; k<N; k++){
            x = k*dx;
            printf("%lf %lf %lf %lf\n", x, ctime-dt, ac[k], bc[k]);
        }

        ac=an;
        bc=bn;
    }
    ac = ac_ptr;
    bc = bc_ptr;

    free(an);
    free(bn);
    free(ac);
    free(bc);

    return 0;
}

int read_coeffs(double *S_arr, double *sigma_arr, int N) {

    FILE *coeff_file;
    coeff_file=fopen("coefficients.txt","r");

    long j;
    for(j=0;j<N;j++) {
        fscanf(coeff_file,"%lf %lf", &S_arr[j], &sigma_arr[j]);
    }
    fclose(coeff_file);

    return 0;
}

void read_input(double *L, int *N, double *D, double *v, double *k_plus, double *k_minus) {
    FILE *infile;
    if(!(infile=fopen("input.txt","r"))) {
        printf("Error opening file\n");
        exit(1);
    }
    if(6!=fscanf(infile,"%lf %d %lf %lf %lf %lf", L, N, D, v, k_plus, k_minus)) {
        printf("Error reading parameters from file\n");
        exit(1);
    }
    fclose(infile);
}

