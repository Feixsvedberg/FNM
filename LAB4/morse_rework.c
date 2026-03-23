#include "numerov.h"
#define tol 1e-16


//Function to calculate morse potential in a point R
double calc_V_morse(double R, Parameters *p)
{
    double V;
    V = p->E_B*gsl_pow_2(1.0-exp(-p->a*(R-p->R_e)));
    return V;
}

//pre-calculate f(R) for the entire domain
//Also return a vector R with radial values
void pre_calc_f(double *f, double *R, Parameters *p)
{
    double Ri;
    double V;

    for(int i = 0; i<p->N; i++)
    {
        Ri = p->R0 + i*p->h;
        R[i] = Ri;
        V = calc_V_morse(Ri, p);
        f[i] = (2*p->mu/gsl_pow_2(p->hbar))*(p->E_n-V);
        //printf("f%d = %.15e", i, f[i]);
    }
}

// function to print the morse potential to get insight
void print_V_morse_to_file(Parameters *p)
{
    double V_morse[p->N];
    double R;
    for(int i = 0; i<p->N; i++)
    {
        R = p->R0 + i*p->h;
        V_morse[i] = calc_V_morse(R, p);
    }
    print_to_file(V_morse, p->N, "output_numerov/V_morse.txt");
}

//Function f to be evaluated by the multiroot finder.
int discrepancy (const gsl_vector *E, void *params, gsl_vector *F)
{
    Parameters *p = (Parameters *)params;
    //Update the energy guess
    p->E_n = gsl_vector_get(E,0);

    int N = p->N;
    double h = p->h;
    //Matching point (does not correspond to an index, just decides where we split the grid)
    int xM = p->xM;
    double f[N];
    double R[N];
    pre_calc_f(f, R, p);
    
    
    //Number of elements from left and right when integrating to one past matching point
    int N_left = xM+2;
    int N_right = N-xM;

    //Create vectors to hold the values of f needed for the two directions
    double f_left[N_left];
    double f_right[N_right];
    for (int i = 0; i<N_left; i++)
    {
        f_left[i] = f[i];
    }
    for (int i = 0; i<N_right; i++)
    {
        f_right[i] = f[xM-1+i];
    }

    //Set initial value and translate them into Y (used by Numerovs)
    double psi_0 = 0;
    double psi_1 = 1e-12;
    double Y0 = psi_to_Y(psi_0, f[0], h);
    double Y1 = psi_to_Y(psi_1, f[1], h);
    double Y1_left = gsl_pow_int(-1,p->n)*Y1;
    double Y_left[N_left];
    double Y_right[N_right];

    //Numerovs method from left and right
    numerov_method(Y_left, N_left, h, f_left, Y0, Y1_left, "fromleft");
    numerov_method(Y_right, N_right, h, f_right, Y0, Y1, "fromright");

    //create vectors containing the psi's needed to calculate discrepancy
    double psi_left[N_left];
    double psi_right[N_right];

    //for(int i = 0; i < 3; i++)
    //{
    //    psi_left[i] = Y_to_psi(Y_left[xM-2+1], f_left[xM-2+i],p->h);
    //    psi_right[i] = Y_to_psi(Y_right[xM-2+i], f_right[xM-2+i], p->h);
    //}

    //Translate the iteration variables back into psi and store them in vectors
    //corresponding to the two directions
    for(int i = 0; i < N_left; i++)
    {
        psi_left[i] = Y_to_psi(Y_left[i], f_left[i],p->h);
    }
    for(int i = 0; i < N_right; i++)
    {
        psi_right[i] = Y_to_psi(Y_right[i], f_right[i],p->h);
    }
    printf("Psi right (denom in scale) %.15e \n", psi_right[1]);
    printf("Psi left %.15e \n", psi_left[N_left-2]);

    //find scale to rescale right solutions to fins match from both direction at xM
    double scale = psi_left[xM]/psi_right[1];
    //rescale right solutions
    for(int i = 0; i<N_right; i++)
    {
        psi_right[i] *= scale;
    }
    

    //calculate discrepancy
    double disc = gsl_pow_2(psi_left[xM-1] - psi_right[0]) + gsl_pow_2(psi_left[xM+1] - psi_right[2]);
    gsl_vector_set(F, 0, disc);

    //debug prints
    printf("Disc: %.15e \n", gsl_vector_get(F,0));
    printf("Energy : %.15e \n", gsl_vector_get(E,0));

    if(disc < tol)
    {
    // Left part 
    for(int i = 0; i < N_left; i++) 
        {
            p->psi[i] = psi_left[i] / R[i];
        }

    // right part 
    int j;
    for(int i = xM+2; i < N; i++) 
    {
        j = i - xM;   
        p->psi[i] = psi_right[j] / R[i];
    } 
            
}
    return GSL_SUCCESS;
}

int main(int argc, char *argv[])
{
    
    //Parameters (p) for particles in morse potential.
    Parameters *p = malloc(sizeof(Parameters));
    if (argc != 2)
    {
        fprintf(stderr, "Enter only one integer energy level\n");
        return 1;
    }
    else
    {
        p->n = atoi(argv[1]);
    }
    
    double m_H = 1.67374*1e-27;
    double m_Cl = 5.88715*1e-26;
    double mu = m_H*m_Cl/(m_Cl+m_H);
    p->mu = mu;
    p->E_B = 7.392*1e-19;
    p->R_e = 1.275*1e-10;
    p->a = 1.812*1e+10;
    p->hbar = GSL_CONST_MKSA_PLANCKS_CONSTANT_HBAR;
    p->omega = sqrt(2*p->E_B*gsl_pow_2(p->a)/p->mu);
    //change here to get initial guess for energy level
    p->E_n = calc_E(p);

    //Grid parameters
    p->R0 = 0.75e-10;
    double RN = 3e-10;
    p->N = 10000;
    p->xM = p->N/4;
    p->h = (RN-p->R0)/(p->N-1);
    p->psi = malloc(p->N*sizeof(double));

    const size_t n_dim = 1;

    gsl_multiroot_function f = {&discrepancy, n_dim, p};

    const gsl_multiroot_fsolver_type *st = gsl_multiroot_fsolver_dnewton;
    gsl_multiroot_fsolver *s = gsl_multiroot_fsolver_alloc(st,n_dim);

    gsl_vector *E = gsl_vector_alloc(n_dim);
    gsl_vector_set(E, 0, p->E_n);

    gsl_multiroot_fsolver_set(s, &f, E);

    int iter = 0;
    int status;
    do
    {
        iter++;
        fprintf(stdout, "iteration %d\n", iter);
        status = gsl_multiroot_fsolver_iterate(s);

        if(status)
        {
            break;
        }
        //Check convergence
        status = gsl_multiroot_test_residual(s->f, tol);
    }
    while (status == GSL_CONTINUE && iter < 1000);

    printf("status = %s \n", gsl_strerror(status));

    //print results
    if(status== GSL_SUCCESS)
    {
        fprintf(stdout, "Eigenenergy = ");
        gsl_vector_fprintf(stdout, s->x, "%.6e");
        normalize(p->psi, p); //does not work here for some reason 
        char psi_filename[100];
        sprintf(psi_filename, "output_numerov/psi_E%d_morse.txt", p->n);
        print_to_file(p->psi, p->N, psi_filename);
    }
    


    //numerov_method(Y_left, N, R0, h, f, Y0, Y1, "fromleft");
    //numerov_method(Y_right, N, R0, h, f, Y0, Y1, "fromright");
    //print_to_file(Y_right, N, "output_numerov/Y_right_check.txt");
    //print_to_file(Y_left, N, "output_numerov/Y_left_check.txt");
    print_V_morse_to_file(p);

    
     
        //print_to_file(R, p->N, "output_numerov/R_morse.txt");
   

    return 0;
}