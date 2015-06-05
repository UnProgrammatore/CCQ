#include "../include/sieve.h"

#include <stdio.h>

unsigned int sieve(
	mpz_t n,
	unsigned int* factor_base,
	unsigned int base_dim,
	pair* solutions,
	unsigned int** exponents,
	mpz_t* As,
	unsigned int poly_val_num,
	unsigned int n_fatt_max
	) {

	mpz_t n_root;
	mpz_t intermed;
	mpz_init(n_root);
	mpz_init(intermed);
	mpz_sqrt(n_root, n);

	unsigned int fact_count = 0;
	unsigned int i, j;
	mpz_t* evaluated_poly;

	init_vector_mpz(& evaluated_poly, poly_val_num);

	// Trovo poly_val_num valori del polinomio (A + s)^2 - n, variando A
	// Temporaneamente tolto per usare calcolo potenzialmente più efficiente
	for(i = 0; i < poly_val_num; ++i) {
		mpz_add_ui(intermed, n_root, i);
		mpz_mul(intermed, intermed, intermed);
		mpz_sub(evaluated_poly[i], intermed, n);

		mpz_add_ui(As[i], n_root, i);
	}

	n_fatt_max += base_dim;

	// Per ogni primo nella base di fattori
	for(i = 0; i < base_dim; ++i) {

		// Provo tutte le possibili fattorizzazioni nella base di fattori
		for(j = solutions[i].sol1; j < poly_val_num; j += factor_base[i]) {

			// Calcolo il polinomio solo se mi serve effettivamente
			/*if(mpz_cmp_ui(evaluated_poly[j], 0) == 0) {
				mpz_add_ui(intermed, n_root, j);
				mpz_mul(intermed, intermed, intermed);
				mpz_sub(evaluated_poly[j], intermed, n);

				mpz_add_ui(As[j], n_root, j);
			}*/
			// Divido e salvo l'esponente va bene
			while(mpz_divisible_ui_p(evaluated_poly[j], factor_base[i])) {
				set_matrix(exponents, j, i, get_matrix(exponents, j, i) + 1); // ++exponents[j][i];
				mpz_divexact_ui(evaluated_poly[j], evaluated_poly[j], factor_base[i]);

			}
			
			if(mpz_cmp_ui(evaluated_poly[j], 1) == 0) {
				++fact_count;
				if(fact_count >= n_fatt_max) {
					printf("L'ultimo valore trovato è stato calcolato con la sol1 di a = %d dividendo per %d\n", j, factor_base[i]);
					i = base_dim;
					j = poly_val_num; // Doppio break
				}
			}
		}

		// Faccio la stessa cosa con entrambe le soluzioni, a meno che non stia usando 2
		if(factor_base[i] != 2 && i != base_dim && j != poly_val_num) {
			for(j = solutions[i].sol2; j < poly_val_num; j += factor_base[i]) {
				
				/*if(mpz_cmp_ui(evaluated_poly[j], 0) == 0) {
					mpz_add_ui(intermed, n_root, j);
					mpz_mul(intermed, intermed, intermed);
					mpz_sub(evaluated_poly[j], intermed, n);

					mpz_add_ui(As[j], n_root, j);
				}*/

				while(mpz_divisible_ui_p(evaluated_poly[j], factor_base[i])) {
					set_matrix(exponents, j, i, get_matrix(exponents, j, i) + 1); // ++exponents[j][i];
					mpz_divexact_ui(evaluated_poly[j], evaluated_poly[j], factor_base[i]);
				}
				if(mpz_cmp_ui(evaluated_poly[j], 1) == 0) {
					++fact_count;
					if(fact_count >= n_fatt_max) {
						printf("L'ultimo valore trovato è stato calcolato con la sol2 di a = %d dividendo per %d\n", j, factor_base[i]);
						i = base_dim;
						j = poly_val_num; // Doppio break
					}
				}		
			}
		}
	}

	mpz_clear(n_root);
	mpz_clear(intermed);

	remove_not_factorized(exponents, evaluated_poly, As, poly_val_num, base_dim);

	finalize_vector_mpz(&evaluated_poly, poly_val_num);

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
	unsigned int k = 0;

	for(i = 0; i < howmany; ++i) {
		if(mpz_cmp_ui(reduced_q_a[i], 1) == 0) {
		  mpz_set(q_a[k], q_a[i]);
		  for(j = 0; j < primes_num; ++j)
		  	set_matrix(exponents, k, j, get_matrix(exponents, i, j));
		  ++k;
		}
	}

	return k;
}
