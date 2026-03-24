#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_sf_trig.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_sf_hermite.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>


typedef struct
{
    double E_B;
    double m;
    double mu;
    double a;
    double R_e;
    double hbar;
    double omega;
    double E_n;
    double angular_const;
    int n;
    int N;
    double h;
    double R0;
    int xM;
    //parameter to store final psi
    double *psi;

}Parameters;


double psi_to_Y(double psi, double f, double h);

double Y_to_psi(double Y, double f, double h);

void numerov_method(double *Y, int N, double h, double *f, double Y0, double Y1, const char *direction);

void normalize(double *psi, Parameters *p);

void print_to_file(double* array, int len, const char *filename);
