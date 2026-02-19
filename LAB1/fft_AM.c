#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_sf_trig.h>
#include "signal_tools.h"



/*
Function that creates an AM-signal

Input:
- N: Number of elements in signal (nr of pairs)
- delta: Time between samples
- f: carrier frequency
- t0: Start time
- u0: Lower modulating amplitude
- u1: Higher modulating amplitude

Output:
A pointer to an array containing the AM signal
*/
gsl_complex_packed_array AM_signal(int N, double delta, double f, double t0, double u0, double u1)
{
   
    gsl_complex_packed_array signal = malloc(2*N*sizeof(double));

    double t = 0.0;
    //Variable to shift between amplitude modulation of 1 and 3. If -1 we enforce 1 and if 1 enforce 3.
    double amp_shift = 1.0; //start at 1 since it will shift at k = 0.
    double amp;
    for (int k = 0; k < N; k++)
    {
        t = t0 + k*delta;
        //If one period of the carrier wave has passed, we shift the modulating amplitude
        if(k%128 == 0)
        {
            //amp_shift works like an on and off button
            amp_shift = amp_shift*(-1);
        }
        if(amp_shift<0)
        {
            amp = u0;
        }
        else
        {
            amp = u1;
        }
        signal[2*k] = amp*gsl_sf_sin(2*M_PI*f*t);
        signal[2*k+1] = 0;
    }
    return signal;
}


int main ()
{
    //Parameters
    int N = 1024; //Number of samples
    double t_start = 0.0;
    double t_end = 1.0;
    double T = t_end-t_start;
    double delta = T/N; //Time between samples
    double freq_c = 1/(128.0*delta); //carrier frequency
    //Modulating amplitudes of scenarios 1 and 2.
    double u0_1 = 1.0;
    double u1_1 = 3.0;
    double u0_2 = 0.5;
    double u1_2 = 3.1225;
    
    // Create signals for three different scenarios and transform them
    gsl_complex_packed_array signal_1 = AM_signal(N, delta, freq_c, t_start, u0_1, u1_1);
    int transform_status_1 = gsl_fft_complex_radix2_forward (signal_1, 1, N);

    gsl_complex_packed_array signal_2 = AM_signal(N, delta, freq_c, t_start, u0_2, u1_2);
    int transform_status_2 = gsl_fft_complex_radix2_forward (signal_2, 1, N);

    gsl_complex_packed_array signal_3 = AM_signal(N, delta, freq_c*2, t_start, u0_2, u1_2);
    int transform_status_3 = gsl_fft_complex_radix2_forward (signal_3, 1, N);


    //reorder FFT signals
    reorder(signal_1, 2*N); reorder(signal_2, 2*N); reorder(signal_3, 2*N);


    // Get frequency spectrums
    double *freq_spectrum_1 = get_freq_spectrum(signal_1, N, delta);
    double *freq_spectrum_2 = get_freq_spectrum(signal_2, N, delta);
    double *freq_spectrum_3 = get_freq_spectrum(signal_3, N, delta);

    //print to files
    print_to_file(freq_spectrum_1, N, "output_fft/fft_freq_spec_1.txt");
    print_to_file(freq_spectrum_2, N, "output_fft/fft_freq_spec_2.txt");
    print_to_file(freq_spectrum_3, N, "output_fft/fft_freq_spec_3.txt");

    //free allocated memory
    free(signal_1); free(signal_2); free(signal_3);
    free(freq_spectrum_1); free(freq_spectrum_2); free(freq_spectrum_3);
    return 0;
}