#ifndef TRIVIAL_FACT_GUARD
#define TRIVIAL_FACT_GUARD

#include <gmp.h>

// Ritorna il numero di fattorizzazioni trovate
unsigned int trivial_fact(
	mpz_t n, // Il numero
	unsigned int* eratostene, // I fattori che verranno provati
	unsigned int erat_dim // La dimensione del vettore eratostene
	);

#endif // TRIVIAL_FACT_GUARD