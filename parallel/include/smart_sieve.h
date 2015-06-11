#ifndef SIEVE_GUARD
#define SIEVE_GUARD

#include "pair.h"
#include "vector.h"
#include "matrix.h"
#include "quadratic_sieve.h"

#include <gmp.h>

unsigned int smart_sieve(
	mpz_t n, // Il numero da fattorizzare
	unsigned int* factor_base, // La base di fattori
	unsigned int base_dim, // La dimensione della base di fattori
	pair* solutions, // Il vettore contenente le soluzioni all'equazione
	unsigned int poly_val_num, // Il numero di valori di A che si vogliono provare nel polinomio
	unsigned int max_fact, // max_fact + base_dim = massime fattorizzazioni
	unsigned int intervals,
	unsigned int startfrom
	);

#endif // SIEVE_GUARD
