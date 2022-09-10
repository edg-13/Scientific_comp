#include <stdio.h>
#include <stdlib.h>

void happy_meter(int size);

int main(){
    happy_meter(1);
    printf("\n");
    happy_meter(2);
    printf("\n");
    happy_meter(3);
    printf("\n");
    happy_meter(4);
    printf("\n");
    happy_meter(5);
    printf("\n");
    happy_meter(6);
    printf("\n");

    return 0;
}

void happy_meter(int size){
    if (size<0){
        printf("Invalid input");
        return;
    }
    else if (size>10){
        printf("Invalid input");
        return;
    }
    else{
        char emoji[] = "I am happy about this assignment";
        printf("%s", emoji);
        int i;
        for(i=0; i<size; i++){
            printf(" :-)");
        }
        return;

    }
}