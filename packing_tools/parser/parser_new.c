#include<stdio.h>

// header format:
// N = 32 ,L = 13.3313044692576223 ,L1= { 13.2755191332318470 , 0.0000000000000000 }  ,L2= { 0.1431421467822201 , 13.3873242220082184 }  ,P = 0.0010000000000000 ,P0= 0.0010000000000000 ,
//

int N;
long double L;
long double L1x;
long double L1y;
long double L2x;
long double L2y;
long double P;
long double P0;

int read_header(FILE * source)
{
    int ctr = 0;
    char j = '0';
    ctr += fscanf(source, " N = %i ,", &N);

    if (ferror(source)) return ctr;
    ctr += fscanf(source, " L = %Lf ,", &L);
    if (ferror(source)) return ctr;
    ctr += fscanf(source, " L1 = { %Lf , %Lf } ,", &L1x, &L1y);
    if (ferror(source)) return ctr;
    ctr += fscanf(source, " L2 = { %Lf , %Lf } ,", &L2x, &L2y );
    if (ferror(source)) return ctr;
    ctr += fscanf(source, " P = %Lf ,", &P);
    if (ferror(source)) return ctr;
    ctr += fscanf(source, " P0 = %Lf ,", &P0);
    return ctr;
}
