/***********************************************************
 * This program demonstrates pointer declarations and
 * various pointer operations. We demonstrate a difference
 * between a scalar variable and an array, both allocated on
 * the stack. 
 * 
 * To compile: gcc -Wall -std=c99 -o w3ws1 w3ex_1.c
 * *************************************************/
#include <stdlib.h>
#include <stdio.h>

int main(void){
    /* Declare two variables and a pointer of the same type */
    long int  var1 = 1, var2 = 5;
    long int* var1ptr;
    long int  i;
    /* Declare a stack array of 10 elements and a pointer */
    double  arr1[10] = {10,11,12,13,14,15,16,17,18,19};
    double  new_elem;
    double* arr1ptr1;
    double* arr1ptr2;
    double* arr1elem;
    
    /* At this point both pointers hold random values that were present
    on the stack locations where these were allocated.
    In order to assign a pointer a valid address we can use '&' operator. */
    var1ptr = &var1;
    
    /* Print variable values and a pointer */
    printf("\n var1: %ld | var2: %ld | var1 memory address: %p\n",var1,var2,var1ptr);
    
    /* We can access the value of the variable that the pointer points to
    using '*' operator. Following statement assigns var2 with the value of var1,
    using de-referencing. */
    var2 = *var1ptr;
    
    /* Print variable values and a pointer */
    printf("\n var1: %ld | var2: %ld | var1 memory address: %p\n",var1,var2,var1ptr);

    /* An array name is a pointer to the start of the memory block occupied
    by this array. These statements are equivalent, but normally we use the first one. */
    arr1ptr1 = arr1;
    arr1ptr2 = &arr1[0];
    /* Print variable values and a pointer */
    printf("arr1ptr1: %p | arr1ptr2: %p\n",arr1ptr1,arr1ptr2);
    
    /* Index can be used to access elements of an array */
    for(i = 0; i < 10; i++){
        var2 = arr1[i] - 5*var1;
        printf("Iteration %ld .... var2 value: %ld\n", i, var2);
    }
    printf("\n\n");
    
    arr1ptr2 = &arr1[9]; i = 0;
    /* But one can also use the increment of the pointer, since C knows the size
    of each element in bits. */
    for(arr1elem = arr1ptr1; arr1elem <= arr1ptr2; arr1elem++){
        new_elem = *arr1elem - 2.5*var1;
        printf("Iteration %ld ....arr1 element %lf ... new_elem: %lf\n", i, *arr1elem, new_elem);
        i++;
    }

    return(0);
}
