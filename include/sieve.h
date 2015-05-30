#ifndef SIEVE_GUARD
#define SIEVE_GUARD

// Aggiungere l'header con le pair

void sieve(
	mpz_t n,
	unsigned int* factor_base,
	unsigned int base_dim,
	pair* solutions,
	unsigned int** exponents,
	unsigned int** exp_mod2,
	unsigned int* first_bit_1,
	unsigned int* count_1,
	unsigned int* valori_a,
	unsigned int* dependencies,
	unsigned int poly_val_num
	);

#endif // SIEVE_GUARD