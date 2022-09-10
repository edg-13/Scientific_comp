#include <stdio.h>
#include <stdlib.h>

void l33t(int number_of_copies_of_the_letter_e);
int i = -1;
int j = 0;

int main(){
    while (i<0) scanf("%d", &i);
    l33t(i);

    return(0);
    
}


void l33t(int n) {
    int j = 0;
    //Need null character at end of string
    char word[n+3];
    word[n+2] = '\0';
    word[0]='l';
    while(j<n){
        word[j+1]='e';
        j++;
    }
    word[n+1] = 't';
    printf(word);
    return;
    
}