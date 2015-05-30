#ifndef SIEVE_GUARD
#define SIEVE_GUARD

// Aggiungere l'header con le pair

void sieve(
	mpz_t n, // Il numero da fattorizzare
	unsigned int* factor_base, // La base di fattori
	unsigned int base_dim, // La dimensione della base di fattori
	pair* solutions, // Il vettore contenente le soluzioni all'equazione
	unsigned int** exponents, // Il vettore di vettori degli esponenti
	unsigned int** exp_mod2, // PROBABILMENTE REMOVIBILE: vettore di vettori degli esponenti modulo 2
	unsigned int* first_bit_1, // PROBABILMENTE REMOVIBILE: vettore delle posizioni dei primi bit a 1
	unsigned int* count_1, // PROBABILMENTE REMOVIBILE: vettore dei numeri di uni nel vettore mod2
	unsigned int* valori_a, // I valori di A
	unsigned int* dependencies, // FORSE REMOVIBILE: Le dipendenze lineari
	unsigned int poly_val_num // Il numero di valori di A che si vogliono provare nel polinomio
	);

#endif // SIEVE_GUARD