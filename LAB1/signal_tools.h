#include <string.h>
#include <gsl/gsl_fft_complex.h>
/*
Function that prints an array of data to a file.

Input:
- array: a pointer to an array of doubles
- len: length of the pointer that array references
- filename: name of file where data is to be printed
*/
void print_to_file(double* array, int len, const char *filename);


/*
Function to reorder the output data. This function expects twice the length. 
This is because the sorting algorithm becomes easier to understand in my opinion. 

Input:
- signal: a pointer to an array containing the ouput data of the fft. 
- N: Number of elements in the array (both complex and real)
*/
void reorder(gsl_complex_packed_array signal, int N);

/*
Function that takes the real part of the gsl_complex_packed_array and 
returns a pointer to an array containing only these values. The function also rescales
the data by delta.

Input:
- signal: a pointer to the gsl_complex_packed_array.
- N: Number of elements in the array (nr of pairs)
- delta: time between samples

Output:
- A pointer to an array containing the real values of the input array.
*/
double* get_real(gsl_complex_packed_array signal, int N, double delta);

/*
Function that takes the complex part of the gsl_complex_packed_array and 
returns a pointer to an array containing only these values. The function also rescales
the data by delta.

Input:
- signal: a pointer to the gsl_complex_packed_array.
- N: Number of elements in the array (nr of pairs)
- delta: time between samples

Output:
- A pointer to an array containing the complex values of the input array.
*/
double* get_complex(gsl_complex_packed_array, int N, double delta);


/*
Function that takes a signal and returns its frequency spectrum. Also rescales the signal by delta.

Input:
- signal: a pointer to an array containing the ouput data of the fft.
- N: Number of elements in the array (nr of pairs)
- delta: time between samples
*/
double* get_freq_spectrum(gsl_complex_packed_array signal, int N, double delta);