/*********************************************************
* This program solves coupled equations which describs
* the diffusion-damped system with constant coefficients.
* A time splitting technique by Marchuk–Strang is used. 
* For two subproblems, A1 and A2 we solve A1 with time length
* tau=(1/2)*dt then solve the other one, A2, for time length
* dt and solve with A1 again for time length tau=(1/2)*dt.
* 
* The read_input function should not be modified, but the
* use of the function in main() may or may not be correct.
*
* To compile: gcc -Wall -Werror -std=c99 -o assign3 assign3.c -lm
*
* List of identified errors:
* 61 - 68 ...... Changed malloc(nx) to (double *)malloc(sizeof(double)*nx)
* 71 - 72 ...... Changed = to == in if expression
* 143 - 146 .... Deleted & symbols from free statements
* 43 ..... Added & in front of variables
* 46 ..... Changed nx to (nx-1)
* 49 ..... Replaced hardcoded dt so that it satisfies stability condition and has an integer no. of steps to finish
* 114 ..... Replaced sfac with -sfac
* 136 ..... Printed ctime-dt instead of ctime
* 100, 118 .....Changed to loop condition to k < nx-1, and k initialised at 1, manually added boundary points
* 83 ..... Changed loop condition to k<nx
* 24 ..... Include math library
* 29 ..... Changed to int main instead of double main
* 96-7 ..... Take out factor of 2 in 2*dt
* 139 ...... Changed %g format specifier to %lf
* 80 ...... Initialise ctime to 0
* 93 ...... Changed condition to ctime < t_F +dt/2
* 62 ...... Initialise uc and vc pointers so the arrays can be freed at end of program
* 127-133 .... Make tmp pointer to un and make uc point to un, and vc to vn
* 130 ..... Reassign uc and vc to original memory locations (not pointing to un and vn)
*--------+--------------------------------------------
*  Line  |     Brief description of a fix
* Number |
*--------+-------------------------------------------
* Example (not a real error):
*  21 ...... Removed space between(void) and {
********************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define PI  3.141592

void read_input(double *D, double *L, int *nx, double *t_F);

int main(void) {
    /******************************/
    /* Declarations of parameters */
    /******************************/
    /* Number of grid points */                  
    int nx;
    /* Length of domain */
    double L;
    /* Equation coefficients */
    double D;
    /* Length of time to run simulation. */
    double t_F;

    /* Read in from file; */
    read_input(&D, &L, &nx, &t_F);

    /* Grid spacing */
    double dx = L/(nx-1);
    double invdx2 = 1.0/(dx*dx);      
    /* Time step */
    double dt = 0.4*(dx*dx)/D;
    long nsteps = t_F/dt;

    nsteps+=1;
    //printf("%ld \n", nsteps);
    dt = t_F/nsteps;
    //printf("%lf", dt*nsteps);

    /************************************************/
    /* Solution Storage at Current / Next time step */
    /************************************************/
    double *uc, *un, *vc, *vn;
    /* Time splitting solutions */
    double *uts1, *uts2, *vts1, *vts2;
    /* Derivative used in finite difference */
    double deriv;

    /*Pointers to allow freeing at end*/
    double *uc_ptr;
    double *vc_ptr;

    /* Allocate memory according to size of nx */
    uc = (double *)malloc(sizeof(double)*nx);
    uc_ptr=uc;
    un = (double *)malloc(sizeof(double)*nx);
    vc = (double *)malloc(sizeof(double)*nx);
    vc_ptr=vc;
    vn = (double *)malloc(sizeof(double)*nx);
    uts1 = (double *)malloc(sizeof(double)*nx);
    uts2 = (double *)malloc(sizeof(double)*nx);
    vts1 = (double *)malloc(sizeof(double)*nx);
    vts2 = (double *)malloc(sizeof(double)*nx);

    /* Check the allocation pointers */
    if (uc==NULL||un==NULL||vc==NULL||vn==NULL||uts1==NULL||
    uts2==NULL||vts1==NULL||vts2==NULL) {
        printf("Memory allocation failed\n");
        return 1;
    }
    
    int k;
    double x;
    /* Current time */
    double ctime = 0;

    /* Initialise arrays */
    for(k = 0; k < nx; k++) {
        x = k*dx;
        uc[k]  = 1.0 + sin(2.0*PI*x/L);
        vc[k]  = 0.0;
        /* Set other arrays to 0 */
        uts1[k] = 0; uts2[k] = 0;
        vts1[k] = 0; vts2[k] = 0;
    }

    /* Loop over timesteps */ 
    while (ctime < t_F + dt/2){    
        /* Rotation factors for time-splitting scheme. */
        double cfac = cos(dt); //TAKE OUT 2
        double sfac = sin(dt);

        /* First substep for diffusion equation, A_1 */ 
        for (k = 1; k < nx-1; k++) {
            x = k*dx;
            /* Diffusion at half time step. */
            deriv = (uc[k-1] + uc[k+1] - 2*uc[k])*invdx2;
            uts1[k] = uc[k] + D * deriv * 0.5*dt;
            deriv = (vc[k-1] + vc[k+1] - 2*vc[k])*invdx2;
            vts1[k] = vc[k] + D * deriv * 0.5*dt;
        }
        deriv = (uc[nx-2] + uc[1] - 2*uc[0])*invdx2;
        uts1[0] = uc[0] + D * deriv * 0.5*dt;
        uts1[nx-1] = uc[nx-1] + D * deriv * 0.5*dt;

        deriv = (vc[nx-2] + vc[1] - 2*vc[0])*invdx2;
        vts1[0] = vc[0] + D * deriv * 0.5*dt;
        vts1[nx-1] = vc[nx-1] + D * deriv * 0.5*dt;

        /* Second substep for decay/growth terms, A_2 */
        for (k = 0; k < nx; k++) {
            x = k*dx;
            /* Apply rotation matrix to u and v, */
            uts2[k] = cfac*uts1[k] + sfac*vts1[k];
            vts2[k] = -sfac*uts1[k] + cfac*vts1[k];
        }

        /* Thirs substep for diffusion terms, A_1 */
        for (k = 1; k < nx-1; k++) {
            x = k*dx;
            deriv = (uts2[k-1] + uts2[k+1] - 2*uts2[k])*invdx2;
            un[k] = uts2[k] + D * deriv * 0.5*dt;
            deriv = (vts2[k-1] + vts2[k+1] - 2*vts2[k])*invdx2;
            vn[k] = vts2[k] + D * deriv * 0.5*dt;
        }
        deriv = (uts2[nx-2] + uts2[1] - 2*uts2[0])*invdx2;
        un[0] = uts2[0] + D * deriv * 0.5*dt;
        un[nx-1] = uts2[nx-1] + D * deriv * 0.5*dt;

        deriv = (vts2[nx-2] + vts2[1] - 2*vts2[0])*invdx2;
        vn[0] = vts2[0] + D * deriv * 0.5*dt;
        vn[nx-1] = vts2[nx-1] + D * deriv * 0.5*dt;
                    
        /* Increment time. */ 
        ctime += dt;  
        for (k = 0; k < nx; k++ ) {
            x = k*dx;
            printf("%lf %lf %lf %lf\n",ctime-dt,x,uc[k],vc[k]);
        }  
        
        /* Copy next values at timestep to u, v arrays. */    
        double *tmp;
        tmp = vn; 
        vc = tmp;
        tmp = un; 
        uc = tmp;
    } 

    /*Reassign vc abd uc to their original memory locations*/
    uc = uc_ptr;
    vc = vc_ptr;

    /* Free allocated memory */  
    free(uc); free(un);
    free(vc); free(vn);
    free(uts1); free(uts2);
    free(vts1); free(vts2);

    return 0;
}

// The lines below don't contain any bugs! Don't modify them
void read_input(double *D, double *L, int *nx, double *t_F) {
    FILE *infile;
    if(!(infile=fopen("input.txt","r"))) {
        printf("Error opening file\n");
        exit(1);
    }
    if(4!=fscanf(infile,"%lf %lf %d %lf",D,L,nx,t_F)) {
        printf("Error reading parameters from file\n");
        exit(1);
    }
    fclose(infile);
}
