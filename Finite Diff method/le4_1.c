/*****************************************************************
* This program uses Euler's method to find numerical solution for
* the 2D trajectory of a projectile in gravitational field with
* drag. We are looking for the trajectory profile only.
*
* Compile: gcc -Wall -g -std=c99 -o le4_1 le4_1.c
* Run: ./le4_1 >> le4_1.dat
******************************************************************/
# include <stdio.h>
# include <stdlib.h>
# include <math.h>

/* Define gravitational acceleration */
#define g     		9.80665
/* Define max time for computation */
#define MAX_TIME    10

/* Function prototype */
void read_input(double *x0, double *y0, double *vx0, double *vy0, double *gamma, double *dt);

int main(void){

     /* Parameters: time step, air drag coefficient */
     double dt, gamma;
     /* Initial conditions */
     double x0, y0, vx0, vy0;
     /* Current/next time step data */
     double vx_c, vy_c, vx_n, vy_n;
     /* Store trajectory in the array */
     double *x, *y;
     /* Number of points in trajectory for a given final time */
     long int nsteps, i;
     /* Read in from file */ 
     read_input(&x0, &y0, &vx0, &vy0, &gamma, &dt);
     nsteps = MAX_TIME / dt;

     /* Allocated position matrices */
     x = (double *)malloc(sizeof(double)*nsteps);
     y = (double *)malloc(sizeof(double)*nsteps);

     /* Validate the pointers to x,y arrays */
     if( (x == NULL) || (y == NULL) ){
          printf("\n!!! Memory allocation failed. Exiting ...!\n");
          return(1);
     }
     
     /* First point on the trajectory */
     x[0]=x0; y[0]=y0;
     /* Current/next time step velocity values */
     vx_c = vx0; vy_c = vy0;
     vx_n = 0; vy_n = 0;
     
     /* Print the initial conditions */
     printf("%lf %lf %lf\n", 0.0, x[0], y[0]);
     /* Start the main calculations. End when MAX_Time is
        reached or when the object hits the ground, at y==0. */
     for(i = 0; i < nsteps-1; i++){
          /* Update the velocities */
          vx_n = vx_c - gamma*vx_c*dt;
          vy_n = vy_c - gamma*vy_c*dt - g*dt;
          /* Update position matrix using current velocity */
          x[i+1] = x[i] + dt*vx_c;
          y[i+1] = y[i] + dt*vy_c;
          /* Copy updated velocity to a current one
          for next iteration of the loop */
          vx_c = vx_n; vy_c = vy_n;
	       /* Print data */
          printf("%lf %lf %lf\n", (i+1)*dt, x[i+1], y[i+1]);
     }

     free(x); free(y);     
     return(0);
}                            

/* This function reads in the initial conditions and the parameters 
   From the file. */ 
void read_input(double *x0, double *y0, double *vx0, double *vy0, double *gamma, double *dt) {
   FILE *infile;
   if(!(infile=fopen("le4_1_input.txt","r"))) {
       printf("Error opening file\n");
       exit(1);
   }
   if(6!=fscanf(infile,"%lf %lf %lf %lf %lf %lf",x0,y0,vx0,vy0,gamma,dt)) {
       printf("Error reading parameters from file\n");
       exit(1);
   }
   fclose(infile);
}
