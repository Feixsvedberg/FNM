#include"signal_tools.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_sf_trig.h>
#include <string.h>



void reorder(gsl_complex_packed_array signal, int N)
{
    double temp;
    for(int k = 0; k<N/2; k++)
    {
        //switch elements k and k + N/2 by using temp variable
        temp = signal[k];
        signal[k] = signal[N/2+k];
        signal[N/2+k] = temp;   
    }
}


void print_to_file(double* array, int len, const char *filename)
{
    //open file for writing
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        perror("failed to open file");
        return;
    }
    for(int i = 0; i < len; i++)
    {
        //print with scientific notation, 15 decimals
        fprintf(fp, "%.15e\n", array[i]);
    }
    fclose(fp);
    printf("data exported to %s\n", filename);

}


double* get_real(gsl_complex_packed_array signal, int N, double delta)
{
    //pointer to allocated memory for array 
    double *real_part = malloc(N*sizeof(double));
    for(int k = 0; k<N; k++)
    {
        real_part[k] = delta*signal[2*k];
    }
    return real_part;
}

double* get_complex(gsl_complex_packed_array signal, int N, double delta)
{
    //pointer to allocated memory for array 
    double *complex_part = malloc(N*sizeof(double));
    for(int k = 0; k<N; k++)
    {
        complex_part[k] = delta*signal[2*k+1]; 
    }
    return complex_part;
}

double* get_freq_spectrum(gsl_complex_packed_array signal, int N, double delta)
{
    double *freq_spectrum = malloc(N*sizeof(double));

    for(int k = 0; k<N; k++)
    {
        //calculate frequency spectrum
        freq_spectrum[k] = (gsl_pow_2(delta*signal[2*k]) + gsl_pow_2(delta*signal[2*k +1]));
    }
    return freq_spectrum;
}

