#ifndef SIEVE_GUARD
#define SIEVE_GUARD

// Aggiungere l'header con le pair

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
	unsigned int* reduced_q_a,
	unsigned int* q_a,
	unsigned int howmany
	);

#endif // SIEVE_GUARD