#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_sf_trig.h>
#include <string.h>
#define dim 3


//FIX HEADER LATER
void frequencies(double A, double B, double m, double *q, double *omega, double *eps);

typedef struct 
{
    const char *name;
    double r;
    double sigma;
    double eps;
    double m;
    const char *quant;
} Material;


//function that calculates A from parameters
double get_A(Material mat)
{
    return(12*mat.eps/gsl_pow_2(mat.r))*(13.0*gsl_pow_int(mat.sigma/mat.r, 12) - 7.0*gsl_pow_int(mat.sigma/mat.r,6));
}

//function that calculates B from parameters
double get_B(Material mat)
{
    return(12*mat.eps/gsl_pow_2(mat.r))*(gsl_pow_int(mat.sigma/mat.r, 6) - gsl_pow_int(mat.sigma/mat.r,12));
}

//function that takes the input material and returns the proper parameters

void set_material_parameters(Material *mat)
{
    double angstrom = 1e-10;
    double eps_scale = 1e-20;
    double mass_scale = 1e-25;

    //Neon
    if(strcmp(mat->name, "Ne")==0)
    {
        mat->sigma = 3.035*angstrom;
        mat->eps = 0.0721*eps_scale;
        mat->r = 3.1562*angstrom;
        mat->m = 0.335092*mass_scale;
    }

    //Argon
    else if(strcmp(mat->name, "Ar")==0)
    {
        mat->sigma = 3.709*angstrom;
        mat->eps = 0.236*eps_scale;
        mat->r = 3.7477*angstrom;
        mat->m = 0.66335*mass_scale;
    }
    
    //Krypton
    else if(strcmp(mat->name, "Kr")==0)
    {
        mat->sigma = 3.966*angstrom;
        mat->eps = 0.325*eps_scale;
        mat->r = 3.9922*angstrom;
        mat->m = 1.3915*mass_scale;
    }
    
    //Xenon
    else if(strcmp(mat->name, "Xe")==0)
    {
        mat->sigma = 4.318*angstrom;
        mat->eps = 0.458*eps_scale;
        mat->r = 4.3346*angstrom;
        mat->m = 2.18017*mass_scale;
    }
    else
    {
        fprintf(stderr, "Please input valid material name: Ne, Ar, Kr or Xe");
        exit(EXIT_FAILURE);
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
    //printf("data exported to %s\n", filename);
}

//Prints the ouput for omega or gamma.
void print_output(double *q, double *quant, double q_N)
{
    for(int i = 0; i < q_N; i++)
    {
        fprintf(stdout, "%f %f %f %9e %9e %9e \n", q[dim*i], q[dim*i+1], q[dim*i+2], quant[dim*i], quant[dim*i+1], quant[dim*i+2]);
    }
    
}

//Function to get a long vector conataining all qs that are evaluated and also fills a c_vector
void get_q_all(double *q_start, double *q_end, double *q_all, double *c_vector, double q_N)
{
    double x_length;
    double y_length;
    double z_length;
    if(q_end != NULL)
    {
        x_length = q_end[0] - q_start[0];
        y_length = q_end[1] - q_start[1];
        z_length = q_end[2] - q_start[2];
        for(int i = 0; i<q_N; i++)
        {
            q_all[dim*i] = q_start[0] + x_length/(q_N-1)*i;
            q_all[dim*i+1] = q_start[1] + y_length/(q_N-1)*i;
            q_all[dim*i+2] = q_start[2] + z_length/(q_N-1)*i;
            c_vector[i] = sqrt(gsl_pow_2(q_all[dim*i]-q_start[0])+gsl_pow_2(q_all[dim*i+1]-q_start[1])+gsl_pow_2(q_all[dim*i+2]-q_start[2]));
        }
    }
    else
    {
        q_all[0] = q_start[0];
        q_all[1] = q_start[1];
        q_all[2] = q_start[2];
    }
}

void calc_freqs(double *q_all, double *freqs, Material mat, double q_N)
{
    double q_temp[3];
    double A = get_A(mat);
    double B = get_B(mat);
    double *omega = calloc(dim, sizeof(double));
    double *eps = NULL;
    for(int i = 0; i<q_N; i++)
    {
        q_temp[0] = q_all[dim*i];
        q_temp[1] = q_all[dim*i+1];
        q_temp[2] = q_all[dim*i+2];
        frequencies(A, B, mat.m, q_temp, omega, eps);
        freqs[dim*i] = omega[0];
        freqs[dim*i+1] = omega[1];
        freqs[dim*i+2] = omega[2];
    }
}

//Calculate summand in eq 18 which is dulong petits law
double dulong_petit(double T, double freq)
{
    double hbar = GSL_CONST_MKSA_PLANCKS_CONSTANT_HBAR;
    double k_B = GSL_CONST_MKSA_BOLTZMANN;
    double val;
    val = k_B*gsl_pow_2(hbar*freq/(k_B*T))*exp(hbar*freq/(k_B*T))/(gsl_pow_2(exp(hbar*freq/(k_B*T))-1));
    return val;
}

int main(int argc, char *argv[])
{

    //Variable status keeps track of if we want to run a program with gamma/omega or cv. It is then
    //used to switch between conditions for the format
    int status = 0;
    if(strcmp(argv[2],"cv")==0)
    {
        status = 1;
    }
    //Initial format checks 
    if(status == 0 && argc != 6 && argc != 9 && argc != 10)
    {
        fprintf(stderr, "Wrong number of arguments. Please provide either 4, 5, 6, 9 or 10 arguments for omega/gamma evaluation\n");
        return 1;
    }
    else if(status == 1 && argc != 4 && argc !=5 && argc != 6)
    {
        fprintf(stderr, "Wrong number of arguments. Please proivde either 4,5 or 6 arguments for cv.");
        return 1;
    }
    //fprintf(stdout, "Valid number of arguments \n");
    //ADD CHECK TO SEE IF RIGHT NUMBER OF INPUT ARGUMENTS

    //Create material
    Material mat; 
    mat.name = argv[1];
    mat.quant = argv[2];


    //ADD CHECK TO SEE INPUT ARGUMENTS: ALSO CREATE FUNCTION TO GET FREQ FOR DIFFERENT Qs

    //Get q-vector from input (case where only one point is specied)

    //Get parameters corresponding to input material
    set_material_parameters(&mat);

    //check if gamma/omega.
    if (status == 0)
    {
        //allocate vector holding starting point
        double *q_start = calloc(dim, sizeof(double));
        double q_N;
        double *q_end = calloc(dim, sizeof(double));

        for(int i = 0; i<dim; i++)
        {
            q_start[i] = atof(argv[3+i]);
        }

        //Check if end point values for q are inputted
        if(argc == 10)
        {
            q_N = atof(argv[9]);
            //Check if valid number of N is inputted.
            if(q_N<=0)
            {
                fprintf(stderr, "Invalid number of points entered. Exiting \n");
                return 1;
            }

            
            for (int i = 0; i < dim; i++)
            {
                q_end[i] = atof(argv[6+i]);
            }
        }
        else if(argc == 9)
        {
            q_N = 11;
            //NU FINNS DENNA I BÅDA: INTE BRA.
            for (int i = 0; i < dim; i++)
            {
                q_end[i] = atof(argv[6+i]);
            }
        }
        // Case where noting more than q_start is inputted.
        else if(argc == 6)
        {
            q_N = 1;
            q_end = NULL;
        }

        double *q_all = calloc(dim*q_N, sizeof(double));
        double *c_vector = calloc(q_N, sizeof(double));

        get_q_all(q_start, q_end, q_all, c_vector, q_N);

        if(strcmp(argv[2], "omega")==0)
        {
            
            //vector to store frequencies in
            double *freqs = calloc(q_N*dim, sizeof(double));
            calc_freqs(q_all, freqs, mat, q_N);

            print_output(q_all, freqs, q_N);
            
            print_to_file(freqs, q_N*dim, "output_phonons/frequencies100.txt");
            print_to_file(c_vector, q_N, "output_phonons/c100.txt"); 
        }

        else if(strcmp(argv[2], "gamma")==0)
        {
            if(q_start[0] == 0 && q_start[1] == 0 && q_start[2] == 0)
            {
                fprintf(stderr, "q = (0,0,0) not defined for gamma");
                return 1;
            }
            double delta = mat.r/1e5;
            //create compressed and expanded material
            Material mat_exp;
            Material mat_comp;
            mat_exp.name = argv[1];
            mat_comp.name = argv[1];
            set_material_parameters(&mat_exp);
            set_material_parameters(&mat_comp);

            //change their nearest neighbor distance r
            mat_exp.r = mat.r + delta;
            mat_comp.r = mat.r - delta;

            double *freqs_exp = calloc(q_N*dim, sizeof(double));
            double *freqs_comp = calloc(q_N*dim, sizeof(double));
            //vector to keep all values of q;
            calc_freqs(q_all, freqs_exp, mat_exp, q_N);
            calc_freqs(q_all, freqs_comp, mat_comp, q_N);

            double *grunesien = calloc(q_N*dim, sizeof(double));

            //Calculate grunesien for all branches of q
            for(int i = 0; i < q_N*dim; i++)
            {
                //printf("Freq exp: %f Freq comp: %f", freqs_exp[i], freqs_comp[i]);
                grunesien[i] = -(log(freqs_exp[i])-log(freqs_comp[i]))/(3*(log(mat_exp.r)-log(mat_comp.r)));
            }
            print_output(q_all, grunesien, q_N);
        }
    }
    //Enters this statement if we are running program for cv
    else
    {
        int T_N;
        double T_end;
        double T_start = atof(argv[3]);

        if (argc == 6)
        {
            T_N = atoi(argv[5]);
            T_end = atof(argv[4]);
            if (T_N <= 0)
            {
                fprintf(stderr, "T_N <= 0, not valid");
                return 1;
            }
        }
        else if(argc == 5)
        {
            T_N = 11;
            T_end = atof(argv[4]);
        }
        else if(argc == 4)
        {
            T_N = 1;
            T_end = T_start;
        }
    

        //qvekt is a file of 48 rows and 4 columns¨
        int rows = 48;
        int cols = 4;
        FILE *fp = fopen("qvekt", "r");
        if (!fp)
        {
            fprintf(stderr, "failed to open qvekt");
            return 1;
        }
        double *qvekt_data = calloc(rows*cols, sizeof(double));

        for(int i = 0; i < rows*cols; i++)
        {
            if(fscanf(fp, "%lf", &qvekt_data[i])!=1)
            {
                fprintf(stderr, "error reading a certain element from file");
                return 1;
            }
        }
        fclose(fp);
        double *q_all = calloc(rows*dim, sizeof(double));
        double *weights = calloc(rows, sizeof(double));

        //put the read data into two separate vectors
        for(int i = 0; i<rows; i++)
        {
            q_all[i*dim] = qvekt_data[i*cols];
            q_all[i*dim+1] = qvekt_data[i*cols+1];
            q_all[i*dim+2] = qvekt_data[i*cols+2];
            weights[i] = qvekt_data[i*cols+3];
        }


        double *freqs = calloc(dim*rows, sizeof(double));
        //calc freqs but with rows as input for q_N since we are calculating for that many qs.
        calc_freqs(q_all, freqs, mat, rows);

        double a = M_SQRT1_2*mat.r;
        double qv_factor = 1.0/8.0*(4.0/1000.0)*gsl_pow_3(1.0/a);

        double T;
        double sum;
        double part_sum;
        double cv_over_V;
        //interval divisor
        for(int i = 0; i<T_N; i++)
        {
            T = T_start + (T_end-T_start)/(T_N-1)*i;
            if (T_N == 1)
            {
                T = T_start;
            }
            sum = 0;
            //q-sum
            for (int q = 0; q<rows; q++)
            {
                part_sum = 0;
                //j-sum
                for(int j = 0; j<dim; j++)
                {
                    part_sum = part_sum + dulong_petit(T, freqs[3*q+j]);
                }
                sum += weights[q]*part_sum;
            }
            cv_over_V = qv_factor*sum;
            fprintf(stdout, "%f %15e \n", T, cv_over_V);
        }
    
    }

    return 0;
    
}