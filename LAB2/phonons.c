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
    printf("data exported to %s\n", filename);
}


void calc_freqs(double *q_start, double *q_end, double q_N, double *freq_mat, double *c_vector, Material mat)
{
    double x_length; double y_length; double z_length;
    if(q_end != NULL)
    {
        x_length = q_end[0] - q_start[0];
        y_length = q_end[1] - q_start[1];
        z_length = q_end[2] - q_start[2];
    }
    //Handle case where q_end is not given
    else
    {
        x_length = 0; 
        y_length = 0; 
        z_length = 0;
    }
    double A = get_A(mat);
    double B = get_B(mat);
    //Below works if it is okay to change the contents every iteration

    double *q = calloc(dim, sizeof(double));
    double *omega = calloc(dim, sizeof(double));
    double *eps = NULL;
    for(int i = 0; i < q_N; i++)
    {
        q[0] = q_start[0] + x_length/(q_N-1)*i;
        q[1] = q_start[1] + y_length/(q_N-1)*i;
        q[2] = q_start[2] + z_length/(q_N-1)*i;
        frequencies(A, B, mat.m, q, omega, eps);
        freq_mat[dim*i] = omega[0];
        freq_mat[dim*i+1] = omega[1];
        freq_mat[dim*i+2] = omega[2];
        c_vector[i] = sqrt(gsl_pow_2(q[0]-q_start[0])+gsl_pow_2(q[1]-q_start[1])+gsl_pow_2(q[2]-q_start[2]));
    }
}


int main(int argc, char *argv[])
{
    //Initial format checks 
    if(argc != 6 && argc != 9 && argc != 10)
    {
        fprintf(stderr, "Wrong number of arguments. Please provide either 6, 9 or 10 arguments\n");
        return 1;
    }
    fprintf(stdout, "Valid number of arguments \n");




    //ADD CHECK TO SEE IF RIGHT NUMBER OF INPUT ARGUMENTS

    //Create material
    Material mat; 
    mat.name = argv[1];


    //ADD CHECK TO SEE INPUT ARGUMENTS: ALSO CREATE FUNCTION TO GET FREQ FOR DIFFERENT Qs

    //Get q-vector from input (case where only one point is specied)

    //Get parameters corresponding to input material
    set_material_parameters(&mat);

    //allocate vector holding starting point
    double *q_start = calloc(dim, sizeof(double));
    double q_N;
    double *q_end = calloc(dim, sizeof(double));
    //
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

        printf("%f \n", q_N);
        
        for (int i = 0; i < dim; i++)
        {
            q_end[i] = atof(argv[6+i]);
            printf("%f \n", q_end[i]);
        }
    }
    else if(argc == 9)
    {
        double q_N = 11;
        //NU FINNS DENNA I BÅDA: INTE BRA.
        for (int i = 0; i < dim; i++)
        {
            q_end[i] = atof(argv[6+i]);
            printf("%f \n", q_end[i]);
        }
    }
    // Case where noting more than q_start is inputted.
    else if(argc == 6)
    {
        q_N = 1;
        q_end = NULL;
    }





    if(strcmp(argv[2], "omega")==0)
    {
        double *c_vector = calloc(q_N, sizeof(double));
        //vector to store frequencies in
        double *freq_save = calloc(q_N*dim, sizeof(double));
        calc_freqs(q_start, q_end, q_N, freq_save, c_vector, mat);
        print_to_file(freq_save, q_N*dim, "output_phonons/frequencies100.txt");
        print_to_file(c_vector, q_N, "output_phonons/c100.txt"); 
    }

    else if(strcmp(argv[2], "gamma"))
    {
        double delta = mat.r/1e3;
        //create compressed and expanded material
        Material mat_exp;
        Material mat_comp;
        set_material_parameters(&mat_exp);
        set_material_parameters(&mat_comp);

        //change their nearest neighbor distance r
        mat_exp.r = mat.r + delta;
        mat_comp.r = mat.r - delta;

        double *freq_save_exp = calloc(q_N*dim, sizeof(double));
        double *c_vector_exp = calloc(q_N, sizeof(double));
        double *freq_save_comp = calloc(q_N*dim, sizeof(double));
        double *c_vector_comp = calloc(q_N, sizeof(double));
        calc_freqs(q_start, q_end, q_N, freq_save_exp, c_vector_exp, mat_exp);
        calc_freqs(q_start, q_end, q_N, freq_save_comp, c_vector_comp, mat_comp);

        double *grunesien = calloc(q_N*dim, sizeof(double));

        for(int i = 0; i < q_N*dim; i++)
        {
            grunesien[i] = -(log(freq_save_exp[i])-log(freq_save_comp[i]))/(3*(log(mat_exp.r)-log(mat_comp.r)));
        }

    }
    return 0;

}