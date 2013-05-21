#include<stdio.h>

// header format:
// N = 32 ,L = 13.3313044692576223 ,L1= { 13.2755191332318470 , 0.0000000000000000 }  ,L2= { 0.1431421467822201 , 13.3873242220082184 }  ,P = 0.0010000000000000 ,P0= 0.0010000000000000 ,
//

struct s_header
{
    int N;
    long double L;
    long double L1x;
    long double L1y;
    long double L2x;
    long double L2y;
    long double P;
    long double P0;
};

int read_header(FILE * source, struct s_header *header)
{
    while(fgetc(source) == '#') while(fgetc(source) != '\n');
    fseek(source, -1, SEEK_CUR);
    int ctr = 0;
    ctr += fscanf(source, " N = %i ,", &header->N);
    if (ferror(source)) return ctr;
    ctr += fscanf(source, " L = %Lf ,", &header->L);
    if (ferror(source)) return ctr;
    ctr += fscanf(source, " L1 = { %Lf , %Lf } ,", &header->L1x, &header->L1y);
    if (ferror(source)) return ctr;
    ctr += fscanf(source, " L2 = { %Lf , %Lf } ,", &header->L2x, &header->L2y);
    if (ferror(source)) return ctr;
    ctr += fscanf(source, " P = %Lf ,", &header->P);
    if (ferror(source)) return ctr;
    ctr += fscanf(source, " P0 = %Lf ,", &header->P0);
    return ctr;
}

int read_particles(FILE * source, long double * storage) {
    fscanf(source, " { ");
    while(fscanf(source, " %Lf ,", storage++));
    fscanf(source, " } ");
}
