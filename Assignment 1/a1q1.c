#include <stdio.h>
#include <stdlib.h>

int isPrime(long k);

int main(){
    printf("%d \n", isPrime(343));
    printf("%d \n", isPrime(1001));
    printf("%d \n", isPrime(123456789));
    printf("%d \n", isPrime(1011289));
    return 0;

}

int isPrime(long k){
    int i;
    for(i=2; i<k; i++){
        if (k%i==0) return 0;
    }
    return 1;
}