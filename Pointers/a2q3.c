#include <stdio.h>
#include <stdlib.h>

//-x c -Wall -Werror -std=c99 -lm

int arr_dot_product(int *arr1, int *arr2, int len);

int main(){
    int a[5] = {-1,2,7,2,5};
    int b[5] = {5,13,-4,2,1};

    printf("Dot product is: %d \n", arr_dot_product(a,b,5));
    return 0;
}

int arr_dot_product(int *arr1, int *arr2, int len){
    int result = 0;
    int i;
    for(i=0; i<len; i++){
        result += arr1[i]*arr2[i];
    }
    return result;
}