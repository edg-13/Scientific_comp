#include <stdio.h>
#include <stdlib.h>

//-x c -Wall -Werror -std=c99 -lm

int  *reverse_order(int *arr, int arr_length);

int main(){
    int a[5] = {4,7,3,4,5};
    int *aptr = a;
    int *test;

    test = reverse_order(aptr, 5);
    printf("%d \n", test[0]);
    printf("%d \n", test[1]);
    printf("%d \n", test[2]);
    printf("%d \n", test[3]);
    printf("%d \n", test[4]);
    
}

int  *reverse_order(int *arr, int arr_length){
    int *tmparr;
    int i;
    tmparr = (int *)(malloc(sizeof(int)*arr_length));
    for(i=arr_length-1; i>=0; i--){
        tmparr[arr_length-1-i] = arr[i];
    }
    
    return tmparr;
}