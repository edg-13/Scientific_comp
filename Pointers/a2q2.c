#include <stdio.h>
#include <stdlib.h>

//-x c -Wall -Werror -std=c99 -lm


int  *cumsum(int *arr, int arr_length);

int main(){
    int a[5] = {4,7,3,4,5};
    int *aptr = a;
    int *test;

    test = cumsum(aptr, 5);
    printf("%d \n", test[0]);
    printf("%d \n", test[1]);
    printf("%d \n", test[2]);
    printf("%d \n", test[3]);
    printf("%d \n", test[4]);
    printf("%d \n", test[5]);
    
}


int  *cumsum(int *arr, int arr_length){
    int *result;
    result = (int *)(malloc(sizeof(int)*(arr_length+1)));
    result[0] = 0;
    int i;
    for(i=0; i<arr_length; i++){
        result[i+1] = result[i] + arr[i];
    }
    return result;
}