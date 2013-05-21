#include <stdlib.h>
#include "parser.c"

int main() {
    FILE * fptr;
    
	struct s_header test;
	
    fptr = fopen("N32~P1e-3~9000.txt", "r");
    
	printf("retval = %i\n", read_header(fptr, &test));
printf("N=%i\n", test.N);
printf("%f\n", test.L1x);
printf("%f\n", test.L1y);
printf("%f\n", test.L2x);
printf("%f\n", test.L2y);
printf("%f\n", test.P);
printf("%f\n", test.P0);
    fclose(fptr);
}