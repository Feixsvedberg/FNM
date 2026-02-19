/***
    
    Program to test the installation of GSL.

    If GSL is installed in its "regular" directory, this program should
    compile and run correctly, using the Makefile accompanying this
    source file.

    make gsl_test
    ./gsl_test

***/

#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_const_mksa.h>

int
main ()
{
  printf ("hbar squared = %e\n", 
	  gsl_pow_2 (GSL_CONST_MKSA_PLANCKS_CONSTANT_HBAR));

  return 0;
}
