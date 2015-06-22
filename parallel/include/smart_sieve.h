#ifndef SIEVE_GUARD
#define SIEVE_GUARD

#include "pair.h"
#include "vector.h"
#include "matrix.h"
#include "quadratic_sieve.h"

#include <gmp.h>
#include <omp.h>

/* */
struct mpz_pair {
  mpz_t sol1;
  mpz_t sol2;
};
typedef struct mpz_pair mpz_pair;

unsigned int smart_sieve(
	mpz_t n, // Il numero da fattorizzare
	unsigned int* factor_base, // La base di fattori
	unsigned int base_dim, // La dimensione della base di fattori
	pair* solutions, // Il vettore contenente le soluzioni all'equazione
	mpz_t begin, // Inizio della sequenza degli A
	unsigned int interval, // A in [begin, begin + interval]
	unsigned int block_size, // Porzione di [startfrom, M] da calcolare alla volta
	unsigned int max_fact // max_fact + base_dim = massime fattorizzazioni
	);

#endif // SIEVE_GUARD
