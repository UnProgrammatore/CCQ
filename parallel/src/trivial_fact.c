#include "../include/trivial_fact.h"


// Ritorna il numero di fattorizzazioni trovate
unsigned int trivial_fact(
	mpz_t n, // Il numero
	unsigned int* eratostene, // I fattori che verranno provati
	unsigned int erat_dim // La dimensione del vettore eratostene
	) {

	unsigned int i;
	for(i = 0; i < erat_dim; ++i) {
		if(mpz_divisible_ui_p(n, eratostene[i]))
			return eratostene[i];
	}

	return 0;

}