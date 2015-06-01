#include "../include/sieve.h"
#include "../include/vector.h"

unsigned int sieve(
	mpz_t n,
	unsigned int* factor_base,
	unsigned int base_dim,
	pair* solutions,
	unsigned int** exponents,
	mpz_t* evaluated_poly_original,
	unsigned int poly_val_num) {

	mpz_t n_root;
	mpz_t intermed;
	mpz_init(n_root);
	mpz_init(intermed);
	mpz_sqrt(n_root, n);

	unsigned int fact_count = 0;
	unsigned int i, j;
	mpz_t* evaluated_poly;
	init_vector_mpz(evaluated_poly, poly_val_num);

	// Trovo poly_val_num valori del polinomio (A + s)^2 - n, variando A
	for(i = 0; i < poly_val_num; ++i) {
		mpz_add_ui(intermed, n_root, i);
		mpz_mul(intermed, intermed, intermed);
		mpz_sub(evaluated_poly[i], intermed, n);
		mpz_set(evaluated_poly_original[i], evaluated_poly[i]);
	}

	// Per ogni primo nella base di fattori
	for(i = 0; i < base_dim; ++i) {

		// Provo tutte le possibili fattorizzazioni nella base di fattori
		for(j = solutions[i].sol1; j < poly_val_num; j += factor_base[i]) {

			// Divido e salvo l'esponente va bene
			while(mpz_divisible_p_ui(evaluated_poly[j], factor_base[i])) {
				set_matrix(exponents, j, i, get_matrix(exponents, j, i) + 1); // ++exponents[j][i];
				mpz_divexact_ui(evaluated_poly[j], evaluated_poly[j], factor_base[i]);
				
			}
			
			if(mpz_cmp_ui(evaluated_poly[j], 1) == 0)
				++fact_count;
		}

		// Faccio la stessa cosa con entrambe le soluzioni, a meno che non stia usando 2
		if(factor_base[i] != 2) {
			for(j = solutions[i].sol2; j < poly_val_num; j += factor_base[i]) {

				while(mpz_divisible_p_ui(evaluated_poly[j], factor_base[i])) {
					set_matrix(exponents, j, i, get_matrix(exponents, j, i) + 1); // ++exponents[j][i];
					mpz_divexact_ui(evaluated_poly[j], evaluated_poly[j], factor_base[i]);
				}
				if(mpz_cmp_ui(evaluated_poly[j], 1) == 0)
					++fact_count;
			
			}
		}
	}

	remove_not_factorized(exponents, evaluated_poly, evaluated_poly_original, poly_val_num, base_dim);

	return fact_count;
}

unsigned int remove_not_factorized(
	unsigned int** exponents,
	mpz_t* reduced_q_a,
	mpz_t* q_a,
	unsigned int howmany,
	unsigned int primes_num
	) {

	unsigned int i;
	unsigned int j;
	unsigned int k;

	for(i = 0; i < howmany; ++i) {
		if(mpz_cmp_ui(reduced_q_a[i], 1) == 0) {
			mpz_set(q_a[k], q_a[i]);
			for(j = 0; j < primes_num; ++j)
				set_matrix(exponents, k, j, get_matrix(exponents, i, j))
			++k;
		}
	}

	return k;
}