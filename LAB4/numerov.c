#include "numerov.h"




//translation from psi to Y
double psi_to_Y(double psi, double f, double h)
{
    double Y = psi + gsl_pow_2(h)/12*f*psi;
    return Y;

}
//translation from Y to psi
double Y_to_psi(double Y, double f, double h)
{
    double psi = Y/(1+gsl_pow_2(h)/12*f);
    return psi;
}

//Numerovs method that takes direction as input and changes how it iterates
void numerov_method(double *Y, int N, double h, double *f, double Y0, double Y1, const char *direction)
{
    //From left
    if (strcmp(direction, "fromleft")==0)
    {
        Y[0] = Y0;
        Y[1] = Y1;
        //have to consider that initial values are already filled, thus starting from 1.
        //Therefore N is shifted so that the last iteration does not try to access an index not in array.
        for(int i = 1; i<N-1; i++)
        {
            Y[i+1] = Y[i]*(2-(gsl_pow_2(h)*f[i])/(1+(gsl_pow_2(h)/12)*f[i]))-Y[i-1];
        }
    }
    //From right
    else if (strcmp(direction, "fromright")==0)
    {
        Y[N-1] = Y0;
        Y[N-2] = Y1;
        for (int i = N-2; i>0; i--)
        {  
            Y[i-1] = Y[i]*(2-(gsl_pow_2(h)*f[i])/(1+(gsl_pow_2(h)/12)*f[i]))-Y[i+1];
        }
    }
    else
    {
        fprintf(stderr, "Enter valid direction: fromleft or fromright \n");
        exit(EXIT_FAILURE);
    }
  
}

//Function that normalizes a wave function psi
void normalize(double *psi, Parameters *p)
{
    double sum = 0;
    int N = p->N;
    double h = p->h;
    for (int i = 1; i<N-1; i++)
    {
        sum = sum + psi[i]*psi[i];
    }
    //normalization constant
    double norm_const;
    norm_const = sqrt(1/(h*(0.5*gsl_pow_2(psi[0]) + sum + 0.5*gsl_pow_2(psi[N-1]))));
    
    for(int i = 0; i<N; i++)
    {
        psi[i] = norm_const*psi[i];
    }
}

//Function that prints array of length len onto file with filename
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