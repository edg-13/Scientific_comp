#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//-x c -Wall -Werror -std=c99 -lm

void en_pot(double *posx, double *posy, double *posz,  long ncharges, double *res);

int main(){
    double a[3] = {1,1,1};
    double b[3] = {-1,1,1};
    double c[3] = {0,0,1};
    double ans;

    en_pot(a, b, c, 3, &ans);
    printf("Answer: %lf \n", ans);

    return 0;
}

void en_pot(double *posx, double *posy, double *posz,  long ncharges, double *res){
    int i;
    int j;
    double P = 0;

    for(i=0; i<ncharges; i++){
        double norm;
        double xcor;
        double ycor;
        double zcor;
        for(j=0; j<ncharges; j++){
            if (i!=j){
                xcor = posx[i] - posx[j];
                ycor = posy[i] - posy[j];
                zcor = posz[i] - posz[j];
                norm = sqrt((xcor*xcor) + (ycor*ycor) + (zcor*zcor));
                P+=0.5/norm;
                printf("P: %lf \n", P);
            }

        }
    }
    
    *res = P;
    printf("Result: %lf \n", *res);

}

