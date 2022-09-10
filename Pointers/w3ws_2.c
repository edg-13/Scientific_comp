/***********************************************************
* Illustrates how to pass pointers to the external function.
* This program uses function swap to arrange two given
* numbers in the ascending order.
*
* gcc -Wall -std=c99 -o w3ws2 w3ws_2.c
****************************************************/

#include <stdio.h>
#include <stdlib.h>

/* Prototype your function here */
void compare(double *pnum1, double *pnum2, int* answer);

int main(void)
{
    int     answer;
    double  num1, num2;
	
    /* Use scanf function to read in two numbers */
    printf("Enter 2 real numbers separated by space: ");
    scanf("%lf %lf", &num1, &num2);

	/* Function call */
	compare(&num1, &num2, &answer);
	
	/* print out result */
	if(answer == 1)
	   printf("Your input in descending order: %lf, %lf\n", num2, num1);
	else
	   printf("Your input in descending order: %lf, %lf\n", num1, num2);
	
    return(0);
}

void compare(double* pnum1, double* pnum2, int* pret)
{
      /* Use dereferencing to read/write to these memory locations */
      *pret = (int)(*pnum1 < *pnum2);
}
