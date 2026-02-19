#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_sf_trig.h>
#include "signal_tools.h"

 /*
 Function to read signal from file

 Input: 
 - filename: name of file to read from
 - len: length of file

 Output:
 pointer to an array contaning the contents of the file
 */

gsl_complex_packed_array read_signal_from_file(const char *filename, int len)
{
    gsl_complex_packed_array signal = malloc(2*len* sizeof(double));
    FILE *fp = fopen(filename, "r");

    //iterate through one line at a time and print contents to array
    for (int i = 0; i<len; i++)
    {
        fscanf(fp, "%lf", &signal[2*i]);
        //set complex part to zero
        signal[2*i+1] = 0;
    }

    fclose(fp);
    return signal;
}
/*
Function that applies a Gaussian filter to the signal

Input:
- signal: pointer to array containing the transformed signal
- N: number of elements in the signal (nr of pairs)
- f_c: carrier frequency
- delta: time between samples

*/
void filter(gsl_complex_packed_array signal, int N, int f_c, double delta)
{
    //filter width
    double width = 150;
    for(int k = 0; k<N; k++)
    {
        //adjust n to take negative frequencies in account
        double n = k - N/2;
        double f = n/(N*delta); //frequency
        signal[2*k] = signal[2*k]*exp(-0.5*gsl_pow_2((fabs(f)-f_c)/width));
        signal[2*k+1] = signal[2*k+1]*exp(-0.5*gsl_pow_2((fabs(f)-f_c)/width));
    }
}

/*
Function that translates an array containing zeros and ones in a bit pattern

Input:
- bit_array: pointer to array of zeros and ones
- len: length of array
*/
double bit_translate(int *bit_array, int len)
{
    double bit_sum = 0;
    for (int i = 0; i<len; i++)
    {
        bit_sum = bit_sum + bit_array[i]*pow(2,len-i-1);
    }
    return bit_sum;
}

/*
Function converts the magnitude of a signal into bits and returns an array of 
ASCII values corresponding to 8 bit characters.

Input:
- signal: pointer to an array of contaning the signal
- N: number of elements in the array.

Output:
An array of of the ASCII values of the 8 bit characters. 
*/
int *bit_reader(double *signal, int N)
{
    int bit_counter = 0;
    double lim = 1.8; //threshold for amplitudes corresponding to 1 (found from investigating signal in matlab).
    
    int pattern_index = 0; //Index of the "bit pattern", i.e the character. 
    int one_read = 0; //Variable to keep track of the amplitude of a 1 has been found
    int bits_per_char = 8;
    int n_chars = 16; // number of characters in the message
    int window = 64; //window to investigate amplitude
    double max = 0; //Tracker of the highest value recorded in a window
    int *bits = malloc(bits_per_char*sizeof(int));
    int *bit_values = malloc(n_chars*sizeof(int)); //Hold ASCII values of the characters to be decoded

    for(int k = 0; k<N; k++)
    {
        if(fabs(signal[k]) > max)
        {
            max = fabs(signal[k]);
        }
        //Check if we are at end of window
        if((k%window == 0 && k>0) || k == N-1)
        {
            //check if our highest recorded value is larger than lim
            if(max>lim)
            {
                one_read = 1;
            }
            bits[bit_counter] = one_read;
            //reset one_read
            one_read = 0;
            bit_counter++;
            //If 8 bits are recorded, we send the array of the ones and zeros to
            //translating function.
            if(bit_counter == 8)
            {
                bit_values[pattern_index] = bit_translate(bits, 8); 
                pattern_index ++;
                //reset bit counter
                bit_counter = 0;
            }
            //reset max
            max = 0; 
        }
    }
    free(bits);
    return(bit_values);
}

/*
Function that takes arrray of ASCII values and prints the corresponding characters

Input:
-bit_values: Array of ASCII values
-len: length of array

*/
void ascii_print(int* bit_values, int len)
{
    for(int i = 0; i < len; i++)
    {
        char character = bit_values[i];
        printf("%c", character);
    }
}

int main (int file, char *filename[])
{
    //parameters
    int N = 8192;
    int f_c = 1024; 
    double delta = 1.0/N;

    //Get signal from file and store it
    gsl_complex_packed_array signal = read_signal_from_file(filename[1], N);

    //transform
    int transform_status = gsl_fft_complex_radix2_forward(signal, 1, N);

    reorder(signal, 2*N);

    double *freq_spectrum = get_freq_spectrum(signal, N, delta);

    print_to_file(freq_spectrum, N, "output_fft/fft_decode_freq_spectrum_1.txt");

    //aplly gaussian filter to transformation
    filter(signal, N, f_c, delta); 

    double *freq_spectrum_filtered = get_freq_spectrum(signal, N, delta);

    print_to_file(freq_spectrum_filtered, N, "output_fft/fft_decode_freq_spectrum_2.txt"); 

    //reorder works in reverse to get intial transformed format back
    reorder(signal, 2*N);

    int reverse_status = gsl_fft_complex_radix2_inverse(signal, 1, N);

    //get real part of the filtered time signal
    double *real_part_filtered = get_real(signal,N,1.0);

    print_to_file(real_part_filtered, N, "output_fft/fft_filtered_time_series.txt");

    //get values that are to be translated into characters
    int *bit_values = bit_reader(real_part_filtered, N);

    //print characters
    ascii_print(bit_values, 16);

    free(signal); free(freq_spectrum); free(freq_spectrum_filtered); 
    free(real_part_filtered); free(bit_values);
}