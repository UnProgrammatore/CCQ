#ifndef SIEVE_GUARD
#define SIEVE_GUARD

// Aggiungere l'header con le pair

#include <gmp.h>

struct pair {
  unsigned int sol1;
  unsigned int sol2;
};

unsigned int sieve(
	mpz_t n, // Il numero da fattorizzare
	unsigned int* factor_base, // La base di fattori
	unsigned int base_dim, // La dimensione della base di fattori
	struct pair* solutions, // Il vettore contenente le soluzioni all'equazione
	unsigned int** exponents, // Il vettore di vettori degli esponenti
	unsigned int poly_val_num // Il numero di valori di A che si vogliono provare nel polinomio
	);

#endif // SIEVE_GUARD
