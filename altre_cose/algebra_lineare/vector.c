#include "vector.h"
#include <gmp.h>

void init_vector(unsigned int** vector, unsigned int n) {
	*vector = (unsigned int*) malloc(sizeof(unsigned int) * n);
}


void finalize_vector(unsigned int** vector) {
	free(*vector);
	*vector = 0;
}

void init_vector_mpz(mpz_t** vector, unsigned int n) {
	*vector = (mpz_t*) malloc(sizeof(mpz_t) * n);
	for(; n > 0; --n)
		mpz_init((*vector)[n - 1]);
}

void finalize_vector_mpz(mpz_t** vector, unsigned int n) {
	for(; n > 0; --n)
		mpz_clear((*vector)[n - 1]);
	free(*vector);
	*vector = 0;
}
