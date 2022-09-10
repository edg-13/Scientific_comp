/***********************************************************
 * This programme solves heat difussion equation using
 * explicit finite difference. Dirichlet boundary conditions
 * are assumed. This is to illustrate the numerical convergence
 * of this simplest explicit method, with varying time step.
 * 
 * Compile: gcc -Wall -g -std=c99 -o heateq heateq.c -lm
 * *********************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int main(void) {
    /* Array to store current and next time step values */
    double pi = 4.0*atan(1.0);
    double *T, *T_n;
    double length = 1.0;
    double deriv, x;
    long nx = 10;
    
    /* Need to include both boundary points at x_0 and x_N,
    so that domain is length 1.0 */
    double dx = length/(nx-1);
    double time_final = 1.0; 
    
    /* Time step set by stability limit dt/dx^2 < 0.5 */
    double dt = 0.5*(dx*dx);
    long nsteps = time_final/dt;
    
    /* In order to finish !exactly! at time_final, increase the
    number of steps by one and decrease the time step slightly:
    stability still satisfied. */
    nsteps += 1;
    dt = time_final/nsteps;

    /* Allocate memory for arrays. Exit if not enough memory. */
    T   = (double *)malloc(sizeof(double)*nx); 
    T_n = (double *)malloc(sizeof(double)*nx); 
    if( (T == NULL) || (T_n == NULL) ){
        printf("!!! Memory allocation failed !!! Exiting...\n");
        return(1);
    }
    
    /* Space and time grid indices respectively */
    long k; long j;
    
    /* Initialisation: delta function at mid-radius
    or sinusoidal initial condition. Can be changed
    by switching the value of init_pulse. Both are
    consistent with the assumed Dirichlet boundary
    conditions T[0] = T[nx-1] = 0.*/
    int init_pulse = 0;
    if (init_pulse){
        for(k = 0; k < nx ;k++){
        T[k] = 0.0;
        }
    T[nx/2] = 1.0;
    }
    else{
        for(k = 0; k < nx; k++){
            x = k*dx;
            T[k] = sin(pi*x);
        }
    }
    /* Print first time step */
    double t0 = 0.0;
    for(k = 0; k < nx ; k++) {
        printf("%g %g %g \n",t0,k*dx,T[k]);
    }
    
    /* Write output parameters to a file to make it easier
     to read in and plot the simulation data. */
     FILE *fp = fopen("heateq_outparams","w");
     fprintf(fp,"%ld %ld %g %g \n", nx, nsteps, dx, dt);
     fclose(fp);

    /* The numerical solution starts here */
    for(j = 0; j < nsteps; j++) {
        /* Loop over interior points, k in [1, nx-1].
        Boundary points always have T[0] = T[nx-1] = 0. */
        for(k = 1; k < (nx-1); k++) {
            deriv = (T[k-1] + T[k+1] - 2*T[k])/(dx*dx);
            T_n[k] = T[k] + deriv * dt;
        }
        T_n[0] = 0;
        T_n[nx-1] = 0;
        /* Copy temperature at next step back into current step */
        for(k = 0; k < nx; k++) {
            T[k] = T_n[k];
            printf("%g %g %g \n",j*dt,k*dx,T[k]);
        }
    }
    
    /* Free allocated memory */
    free(T);
    free(T_n);
    
    return 0;
}
