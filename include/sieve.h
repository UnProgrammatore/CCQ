#ifndef SIEVE_GUARD
#define SIEVE_GUARD

#include <gmp.h>
#include "pair.h"

unsigned int sieve(
	mpz_t n, // Il numero da fattorizzare
	unsigned int* factor_base, // La base di fattori
	unsigned int base_dim, // La dimensione della base di fattori
	pair* solutions, // Il vettore contenente le soluzioni all'equazione
	unsigned int** exponents, // Il vettore di vettori degli esponenti
	mpz_t* evaluated_poly_original, // Vettore in cui verranno salvati i Q(A), inizializzare con l'apposita funzione
	unsigned int poly_val_num // Il numero di valori di A che si vogliono provare nel polinomio
	);

unsigned int remove_not_factorized(
	unsigned int** exponents,
	mpz_t* reduced_q_a,
	mpz_t* q_a,
	unsigned int howmany,
	unsigned int primes_num
	);

#endif // SIEVE_GUARD
