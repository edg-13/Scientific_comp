/**********************************************
* This program demonstrates various pointer
* operations. The compiler should warn you
* about the incorrect pointer use. Try to
* fix this and compare the the bahviour on
* nenneke and your own computer.
*
* Compile: gcc -Wall -g -o w2ex3 w2ex3.c
***********************************************/
#include<stdio.h>
#include<stdlib.h>

int main(void) {

  int *kvec;
  int *jvec;

  int j=5;

  printf("Value of kvec (pointer address) before assignment: %p\n",kvec); 
  printf("Value stored at address pointed to by kvec: %d\n",*kvec);
  kvec = &j;
  printf("Value of kvec (pointer address) after assignment: %p\n",kvec); 
  printf("Value stored at address pointed to by kvec: %d\n",*kvec);

  kvec = (int *)malloc(sizeof(int)*5);

  kvec[3] = 42;
  printf("Value of kvec[3]: %d\n",kvec[3]); 

  jvec = &(kvec[1]);
  printf("Value of jvec[2]: %d ... jvec = &(kvec[1]) \n",jvec[2]); 
  jvec = kvec+1;
  printf("Value of jvec[2]: %d ... jvec = kvec+1 \n",jvec[2]);   


}