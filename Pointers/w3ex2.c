/**********************************************
* This program allocates memory for the aray
* of 5 elements and performs some comparisons
* on pointers.
* 
* Compile: gcc -Wall -g -o w2ex2 w2ex2.c
***********************************************/
#include<stdio.h>
#include<stdlib.h>

int main(void) {

  int *jvec;
  int *kvec;

  /* Explicit casting of the pointer returned from malloc
     is a good programming practice. */
  jvec = (int *)malloc(sizeof(int)*5);
  kvec = (int *)malloc(sizeof(int)*5);

  /* The pointer may have a NULL value, at which point we may
     not want to proceed. */
  if( (jvec == NULL) || (kvec == NULL)){
    printf("Memory allocation failed. Exiting...\n");
    return(-1);
  }
  else{
    printf("Value of jvec (pointer address) after assignment: %p\n",jvec); 
    printf("Value of kvec (pointer address) after assignment: %p\n",kvec);
    printf("Value of NULL (pointer address): %p\n",(int *) NULL); 
    printf("Note on some systems malloc(0) returns NULL!\n");
    
    printf("kvec compared to NULL: (kvec==NULL) %d \n",(kvec==NULL)); 
    printf("jvec compared to kvec: (kvec==jvec) %d \n",(kvec==jvec));
    
    free(jvec);
    free(kvec);
    return(0);
  }  
}