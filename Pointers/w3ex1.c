/***************************************************
* This program demonstrates variables' scope in C.
* Two if(1) statements are just a way of marking
* two different code blocks (enclosed in{}).
* The variables are only visible within this block.
* Pointers allows us to access variables regardless
* their scope.
* If you run this example on Windows 10 machine it
* may take few seconds before the program finishes
* as it gets very confused.
*
* Compile: gcc -Wall -g -o w3ex1 w3ex1.c
***********************************************/
#include<stdio.h>
#include<stdlib.h>

/* This is how NOT to return a pointer to an array */ 
int *return_an_array_incorrectly() {
  int k[0];
  return k;
}

/* There are no comments in the code below, since the
   messages printed to the terminal window explain what
   should be occuring. */
int main(void) {

  int *jvec;
  int *kvec;

  if(1) {
    /* Array allocated on the stack */
    int cake_numbers[5] = {1,2,3,4,5};
    kvec = cake_numbers;
    printf("address of array cake_numbers:                    %p\n",cake_numbers);
    printf("Value of kvec (pointer address) after assignment: %p\n",kvec); 
  }
  
  if(1) {
    int bank_accounts[5] = {11, 12, 13, 14, 15};
    printf("address of array bank_accounts:                   %p\n",bank_accounts);
    printf("Value at bank_accounts[4]: %d\n",bank_accounts[4]);
    printf("Now, assign a value to kvec, which is pointing to an out of scope variable.\n"); 
    bank_accounts[4] = kvec[4];
    printf("Consequences are UNDEFINED!\n"); 
    kvec[4] = 44;
    printf("Value at bank_accounts[4]: %d\n",bank_accounts[4]);
  }

  jvec = return_an_array_incorrectly();
  jvec[4] = 10;
  
  return 0;
}