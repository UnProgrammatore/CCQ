#ifndef QUADRATIC_SIEVE_H
#define QUADRATIC_SIEVE_H 1

#include "../include/linear_algebra.h"
#include "../include/vector.h"
#include "../include/sieve.h"
#include "../include/matrix.h"
#include "../include/base_fattori.h"
#include "../include/trivial_fact.h"

#include <omp.h>
#include <gmp.h>

enum qs_error_codes{
  OK,
  SOLO_FATTORIZZAZIONI_BANALI,
  
};

unsigned long quadratic_sieve(mpz_t N, 
			      unsigned int n, 
			      unsigned int poly_val_num,
			      mpz_t m);

#endif // QUADRATIC_SIEVE_H
