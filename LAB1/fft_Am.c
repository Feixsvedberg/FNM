#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_sf_trig.h>



//Assemble signal
gsl_complex_packed_array AM_signal(int N, double delta, double f, double t0)
{
   
    gsl_complex_packed_array signal = malloc(2*N*sizeof(double));

    double t = 0.0;
    //Variable to shift between amplitude modulation of 1 and 3. If -1 we enforce 1 and if 1 enforce 3.
    double amp_shift = 1.0; //start at 1 since it will shift at k = 0.
    double amp_plus = 0;
    for (int k = 0; k < N; k++)
    {
        t = t0 + k*delta;
        if(k%128 == 0)
        {
            amp_shift = amp_shift*(-1);
        }
        if(amp_shift<0)
        {
            amp_plus = 0;
        }
        else
        {
            amp_plus = 2;
        }
        signal[2*k] = (1+amp_plus)*gsl_sf_sin(2*M_PI*f*t);
        signal[2*k+1] = 0;
    }
    return signal;
}

//Print array to file
void print_to_file(gsl_complex_packed_array array, int len)
{
    const char *filename = "output_fft/fft_AM_data.txt";
    // Open file for writing

    FILE *fp = fopen(filename, "w");
    if (!fp) {
        perror("Failed to open file");
        return;
    }
    for(int i = 0; i < len; i++)
    {
        fprintf(fp, "%f\n", array[i]);
    }
    fclose(fp);
    printf("Data exported to %s\n", filename);
}

int main ()
{
    //Parameters
    int N = 1024; //Number of samples
    double t_start = 0.0;
    double t_end = 1.0;
    double T = t_end-t_start;
    double delta = T/N;
    double freq_c = 1/(128.0*delta);
    
    /****************** */

    gsl_complex_packed_array signal = AM_signal(N, delta, freq_c, t_start);

    //0 if success
    int transform_status = gsl_fft_complex_radix2_forward ((gsl_complex_packed_array)signal, 1, N);

    printf("%d \n", transform_status);

    print_to_file(signal,2*N);
  return 0;
}