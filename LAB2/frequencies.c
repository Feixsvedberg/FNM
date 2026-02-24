#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>


double get_diag_value(double A, double B, double m, double *q, int i)
{
    //adjust indices for permutation
    double A_B_terms; double B_terms;
    int x1 = i; int x2 = i+1; int x3 = i+2;
    if(i == 1) {x3 = 0;}
    else if(i == 2) {x2 = 0; x3 = 1;}
    //Gather terms
    A_B_terms = 8-4*cos(q[x1]*M_PI)*cos(q[x2]*M_PI)-4*cos(q[x1]*M_PI)*cos(q[x3]*M_PI);
    B_terms = 4-4*cos(q[x2]*M_PI)*cos(q[x3]*M_PI);
    return((A+B)*A_B_terms/(2*m) + B*B_terms/m);
}

//Function to calculate the frequencies omega of a given material for a given vector q. 
void frequencies(double A, double B, double m, double *q, double *omega, double *eps)
{
    int dim = 3;
    double diag_val;
    double val;
    
    //Allocate a gsl matrix with zeros and return a pointer to to it
    gsl_matrix *D = gsl_matrix_calloc(dim,dim);
    for(int i = 0; i<dim; i++)
    {
        //set digonal elements
        diag_val = get_diag_value(A,B,m,q,i);
        gsl_matrix_set(D,i,i,diag_val);
        //construct for loop so that it only loops over upper off diagonal matrix, then using symmetry
        for(int j = i+1; j<dim; j++)
        {
                val = (A-B)/(2*m)*(4*sin(q[i]*M_PI)*sin(q[j]*M_PI));
                //set off diagonal elements using symmetry
                gsl_matrix_set(D,i,j,val);
                gsl_matrix_set(D,j,i,val);
        }
    }

    //Obtain eigenvalues, i.e frequencies squared. 
    gsl_eigen_symmv_workspace *work = gsl_eigen_symmv_alloc(dim);
    gsl_vector *eval = gsl_vector_alloc (dim);
    gsl_matrix *evec = gsl_matrix_alloc(dim, dim);
    //calculate eiegenvalues
    gsl_eigen_symmv(D, eval, evec, work);
    //sort vectors
    gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_ASC);


    for(int i = 0; i<dim; i++)
    {

        //write the eigenvalues onto the frequency vector
        double eig_val = gsl_vector_get(eval,i);
        //check if eigenvalues negative and set to zero if that is the case. 
        if(eig_val < 0)
        {
            eig_val = 0;
        }
        omega[i] = sqrt(eig_val);
        //Check if we are interested in eigennvectors
        if (eps != NULL)
        {
            for (int j = 0; j < dim; j++)
            {
                eps[i*dim+j] = gsl_matrix_get(evec,j,i);
            }
        }
    }
    gsl_eigen_symmv_free(work);
    gsl_vector_free(eval);
    gsl_matrix_free(evec);
    gsl_matrix_free(D);

}