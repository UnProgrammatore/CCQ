#ifndef VECTOR_GUARD
#define VECTOR_GUARD
#include <stdlib.h>
#include <gmp.h>


void init_vector(unsigned int** vector, unsigned int n);
void finalize_vector(unsigned int**);

void init_vector_mpz(mpz_t** vector, unsigned int n);
void finalize_vector_mpz(mpz_t** vector, unsigned int n);

#endif // VECTOR_QUARD