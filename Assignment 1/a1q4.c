#include <stdio.h>
#include <stdlib.h>

int num_conseq_digits(long k);

int main(){
    printf("%d", num_conseq_digits(13));
    return 0;
}

int num_conseq_digits(long k){
    long kdup = k;
    int n = 0;
    while(k!=0){
        n++;
        k/=10;
    }
    int digits[n];
    int j;
    for(j=0; j<n; j++){
        digits[j] = kdup%10;
        kdup/=10;
    }
    int count=1;
    int maxcount=1;
    int p;
    for(p=0; p<n-1; p++){
        if(digits[p]==digits[p+1]){
            count++;
            if(count>maxcount){
                maxcount=count;
            }
        }
        else count=1;
    }
    return maxcount;
}