#ifndef SIEVE_GUARD
#define SIEVE_GUARD

#include "pair.h"
#include "vector.h"
#include "matrix.h"

#include <gmp.h>

unsigned int sieve(
	mpz_t n, // Il numero da fattorizzare
	unsigned int* factor_base, // La base di fattori
	unsigned int base_dim, // La dimensione della base di fattori
	pair* solutions, // Il vettore contenente le soluzioni all'equazione
	unsigned int** exponents, // Il vettore di vettori degli esponenti
	mpz_t* As, // Vettore in cui verranno salvati i valori di (A + s)
	unsigned int poly_val_num, // Il numero di valori di A che si vogliono provare nel polinomio
	unsigned int max_fact // max_fact + base_dim = massime fattorizzazioni
	);

unsigned int remove_not_factorized(
	unsigned int** exponents,
	mpz_t* reduced_q_a,
	mpz_t* As,
	unsigned int howmany,
	unsigned int primes_num
	);

#endif // SIEVE_GUARD
