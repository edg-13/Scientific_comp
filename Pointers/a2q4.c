#include <stdio.h>
#include <stdlib.h>

//-x c -Wall -Werror -std=c99 -lm

void varr_swap(char *arr, int a, int b);

int main(){
    char a[5] = {'a','b','c','d','e'};
    varr_swap(a, 0, 4);

    printf("%c \n", a[0]);
    printf("%c \n", a[1]);
    printf("%c \n", a[2]);
    printf("%c \n", a[3]);
    printf("%c \n", a[4]);
    
    return 0;
}

void varr_swap(char *arr, int a, int b){
    char tmp[2];
    tmp[0] = arr[a];
    tmp[1] = arr[b];

    arr[a] = tmp[1];
    arr[b] = tmp[0];
}