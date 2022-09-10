/***********************************************
 * This program solves decay system ODE
 * using explicit finite difference.
 * Time step dt should be varied to see
 * how accurate of a solution we can get.
 * 
 * To compile: gcc -Wall -g -o le4_2 le4_2.c
 * Run (on nenneke): ./le4_2 >> le4_2out.dat
 * This line generates an output file le4_2out.dat
 * which can be analysed (plotted and compared with
 * the analytic solution).
 * **************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int main(void) {
  double time_final = 10.0;
  /* Set decay constant */ 
  double gamma = -4.0;
  /* Set time step */
  double dt = 0.010;
  long nsteps = time_final/dt;
  long i;
 
  /* Initial condition */
  double y = 1.0;
  
  /* The loop gives solution at each time step */
  for(i=0;i<nsteps;i++) {
    printf("%g %g\n",(i*dt),y);
    y = y + gamma*dt*y;     
  }
  printf("%g %g \n",(i*dt),y);

  return(0);
}
