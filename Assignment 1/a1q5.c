#include <stdio.h>
#include <stdlib.h>

long int a;
long int b;

int main(){
    if (scanf("%ld %ld", &a, &b)!=2) {
        printf("Not all field were assigned \n");
    }
    else {
        long int result = abs(a-b);
        printf("%ld", result);
    }
    return 0;
}