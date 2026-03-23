#include "numerov.h"

//Calculate potential
double calc_V(double x, Parameters p)
{
    double V = p.m*gsl_pow_2(p.omega*x)/2.0;
    return V;
}



//pre-calculate f(x)
void pre_calc_f(double *f, double x0, Parameters p)
{
    double h = p.h;
    int N = p.N;
    double xi;
    double V;

    for(int i = 0; i<N; i++)
    {
        xi = x0 + i*h;
        V = calc_V(xi, p);
        f[i] = (2*p.m/gsl_pow_2(p.hbar))*(p.E_n-V);
    }
}

void calc_psi_diff(double *diff, double *psi, double *psi_analytical, int N)
{
    for (int i = 0; i < N; i++)
    {
        diff[i] = psi[i] - psi_analytical[i];
    }

}




void calc_analytical_psi(double *psi_analytical, double x0, Parameters p)
{
    int N = p.N;
    double h = p.h;
    double x; 
    double H_n;
    int n = p.n;
    for(int i = 0; i < N; i++)
    {
        x = x0 + h*i;
        H_n = gsl_sf_hermite(n, x);
        psi_analytical[i] = 1.0/(sqrt(M_SQRTPI)*sqrt(gsl_pow_int(2,n)*gsl_sf_gamma(n+1)))*H_n*exp(-gsl_pow_2(x)/2);
    }
}

int main()
{
    for(int n = 0; n<5; n++)
    {

        //Parameters for a picle (including hbar)
        Parameters p;
        p.n = n;
        p.m = 1.0;
        p.hbar = 1.0;
        p.omega = 1.0;
        p.E_n = calc_E(&p);
        p.N = 1000;
        //N is the size of the two sides
        int N_half = p.N/2;
        double x0 = -8.0;
        double xN = 8.0;
        p.h = (xN-x0)/(p.N-1);
        double f[p.N];
        pre_calc_f(f, x0, p);
        double f_left[N_half];
        double f_right[N_half];
        split_f(f, f_left, f_right, N_half);

        double psi_0 = 0;
        double psi_1 = 1e-12;
        double Y0 = psi_to_Y(psi_0, f_left[0], p.h);
        double Y1 = psi_to_Y(psi_1, f_left[1], p.h);
        //adjust for the odd functions. Symmetry but both negative and positive are accepted solutions
        double Y1_left = gsl_pow_int(-1,p.n)*Y1;

        double Y_left[N_half];
        double Y_right[N_half];

        numerov_method(Y_left, N_half, p.h, f_left, Y0, Y1_left, "fromleft");
        numerov_method(Y_right, N_half, p.h, f_right, Y0, Y1, "fromright");

        double psi[p.N];

        for(int i = 0; i < N_half; i++)
        {
            //printf("Y%d = %5e ", i, Y_left[i]);
            //printf("Y%d = %5e \n", N/2-1-i, Y_right[N/2-1-i]);
            //printf("fleft%d = %5e \n", i, f[i]);
            //printf("fright%d = %5e \n", N-1-i, f_right[N-1-i]);
            psi[i] = Y_to_psi(Y_left[i], f_left[i], p.h);
            psi[N_half+i] = Y_to_psi(Y_right[i], f_right[i], p.h);
        }
        normalize(psi, &p);

        double psi_analytical[p.N];

        calc_analytical_psi(psi_analytical, x0, p);

        double psi_diff[p.N];

        calc_psi_diff(psi_diff, psi, psi_analytical, p.N);

        char psi_filename[100];
        sprintf(psi_filename, "output_numerov/psi_E%d.txt", p.n);
        char psi_anal_filename[100];
        sprintf(psi_anal_filename, "output_numerov/psi_E%d_exact.txt", p.n);
        char diff_filename[100];
        sprintf(diff_filename, "output_numerov/psi_E%d_diff.txt", p.n);

        //print_to_file(Y_left, N, "output_numerov/Y_left.txt");
        //print_to_file(Y_right, N, "output_numerov/Y_right.txt");
        print_to_file(psi, p.N, psi_filename);
        print_to_file(psi_analytical, p.N, psi_anal_filename);
        print_to_file(psi_diff, p.N, diff_filename);
    }
    return 0;
}