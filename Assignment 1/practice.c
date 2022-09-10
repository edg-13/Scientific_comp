#include <stdio.h>
#include <stdlib.h>

long int A;
long int B;
long long int power(long int a, long int b);

int main(void){

    //printf("Enter base followed by exponent \n");

    if (scanf("%ld %ld", &A, &B)!=2) {
        printf("Not all field were assigned \n");
    }
    else {
        if(A==1){
            printf("1");
            return 0;
        }
        if(A==-1){
            if(B%2==0) printf("1");
            else printf("-1");
            return 0;
        }
        long long int result = power(A,B);
        printf("%lld", result);

    return(0);

    }
    
}
long long int power(long int a, long int b) {
    long long int j = 1;
    while(b>0){
        j=j*a;
        b--;
    }
    return(j);
}