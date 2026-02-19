#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_sf_trig.h>
#include "signal_tools.h"

//Gaussian function
double gaussian(double x)
{
    double sigma = 1.0/64.0;
    double val = 1.0/sqrt(M_PI*gsl_pow_2(sigma))*exp(-gsl_pow_2(x/sigma));
    return val;
}

/*
Function to assemble gaussian signal

Input:
- N: Number elements in the signal
- delta: time between samples
- t0: Start time

Output: 
A pointer to a an array containing the signal
*/
gsl_complex_packed_array gauss_signal(int N, double delta, double t0)
{
     //2N because of complex entries
    gsl_complex_packed_array signal = malloc(2*N*sizeof(double));
    //printf("%d \n", N);
    for (int k = 0; k < N; k++)
    {
        double t = t0 + k*delta;
        signal[2*k] = gaussian(t);
        signal[2*k+1] = 0; //complex entries zero
    }
    return signal; 
}


/*
Function to correct for time shift.

Input:
- signal: a pointer to an array containing the ouput data of the fft.
- N: number of elements in the array (nr of pairs)

*/
void time_correction(gsl_complex_packed_array signal, int N)
{
    for(int k=0; k<N; k++)
    {
        //adjust n to take negative frequencies in account
        double n = k - N/2;
        //No sin since pi-integer arg. 
        signal[2*k] = signal[2*k]*cos(M_PI*n);
        signal[2*k+1] = signal[2*k+1]*cos(M_PI*n);
    }
}

/*
Function to get analytical transformation

Input:
- N: Number of real elements in FFT output
- delta: Time between samples

Ouput:
- A pointer to an array of doubles containing the analytical tranformation of a gaussian signal
*/
double *analytical_trans(int N, double delta)
{
    double *arr = malloc(N*sizeof(double));
    double sigma = 1.0/64.0;
    for(int k = 0; k<N; k++)
    {
        //adjust n to take negative frequencies in account
        double n = k - N/2;
        double f = n/(N*delta);//frequencies
        arr[k] = exp(-0.25*gsl_pow_2(2*M_PI*f*sigma));
    }
    return arr;
}



int main ()
{
    //Parameters
    int N = 2048; //Number of samples
    double t_start = -1.0/3.0;
    double t_end = 1.0/3.0;
    double T = t_end-t_start;
    double delta = T/N; // Time between samples


    //Get gaussian signal
    gsl_complex_packed_array signal = gauss_signal(N, delta, t_start);

    //0 if success
    int transform_status = gsl_fft_complex_radix2_forward (signal, 1, N);

    //printf("%d \n", transform_status);

    //Reorder fft ouput so that it is in order of increasing frequencies
    reorder(signal,2*N);

    //Perform time correction since our starting time is not 0
    time_correction(signal,N);

    //get real and complex parts
    double *real_part = get_real(signal,N, delta);
    double *complex_part = get_complex(signal,N, delta);

    //Print real and complex parts to file
    print_to_file(real_part, N, "output_fft/fft_gauss_real.txt");
    print_to_file(complex_part, N, "output_fft/fft_gauss_complex.txt");

    //get analytical solution
    double *exact_trans = analytical_trans(N, delta);

    print_to_file(exact_trans, N, "output_fft/fft_exact_trans.txt");

    //free allocated memory
    free(real_part);
    free(complex_part);
    free(exact_trans);
    free(signal);
  return 0;
}