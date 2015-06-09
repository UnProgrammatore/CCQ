#include "../include/sieve.h"
#include "../include/linear_algebra.h"
#include "../include/matrix.h"

unsigned int sieve(
	mpz_t n,
	unsigned int* factor_base,
	unsigned int base_dim,
	pair* solutions,
	unsigned int** exponents,
	mpz_t* As,
	unsigned int poly_val_num,
	unsigned int max_fact
	) {

	mpz_t n_root;
	mpz_t intermed;
	mpz_init(n_root);
	mpz_init(intermed);
	mpz_sqrt(n_root, n);

	unsigned char go_on = 1;

	unsigned int fact_count = 0;
	unsigned int i, j, k;
	mpz_t* evaluated_poly;

	init_vector_mpz(&evaluated_poly, poly_val_num);

	word** is_used;
	init_matrix_l(&is_used, 1, (poly_val_num / N_BITS) + 1);
	for(i = 0; i < ((poly_val_num / N_BITS) + 1); ++i) {
		set_matrix_l(is_used, 0, i, 0);
	}

	max_fact += base_dim;

	// Trovo poly_val_num valori del polinomio (A + s)^2 - n, variando A
	for(i = 0; i < poly_val_num; ++i) {
		mpz_add_ui(intermed, n_root, i);
		mpz_mul(intermed, intermed, intermed);
		mpz_sub(evaluated_poly[i], intermed, n);
		
		mpz_add_ui(As[i], n_root, i);
	}

	// Per ogni primo nella base di fattori
	for(i = 0; i < base_dim && go_on; ++i) {

		// Provo tutte le possibili fattorizzazioni nella base di fattori
		for(j = solutions[i].sol1; j < poly_val_num && go_on; j += factor_base[i]) {

			// Divido e salvo l'esponente va bene
			while(mpz_divisible_ui_p(evaluated_poly[j], factor_base[i])) {
				
				// Se non sono mai stati usati gli esponenti
				if(get_k_i(is_used, 0, j) == 0) {
					for(k = 0; k < base_dim; ++k)
						set_matrix(exponents, j, k, 0);
					set_k_i(is_used, 0, j, 1);
				}
				
				set_matrix(exponents, j, i, get_matrix(exponents, j, i) + 1); // ++exponents[j][i];
				mpz_divexact_ui(evaluated_poly[j], evaluated_poly[j], factor_base[i]);	
			}
			
			if(mpz_cmp_ui(evaluated_poly[j], 1) == 0) {
				++fact_count;
				if(fact_count >= max_fact) {
					go_on = 0;
				}
			}
		}

		// Faccio la stessa cosa con entrambe le soluzioni, a meno che non stia usando 2
		if(factor_base[i] != 2) {
			for(j = solutions[i].sol2; j < poly_val_num && go_on; j += factor_base[i]) {

				while(mpz_divisible_ui_p(evaluated_poly[j], factor_base[i])) {
					
					// Se non sono mai stati usati gli esponenti
					if(get_k_i(is_used, 0, j) == 0) {
						for(k = 0; k < base_dim; ++k)
							set_matrix(exponents, j, k, 0);
						set_k_i(is_used, 0, j, 1);
					}

					set_matrix(exponents, j, i, get_matrix(exponents, j, i) + 1); // ++exponents[j][i];
					mpz_divexact_ui(evaluated_poly[j], evaluated_poly[j], factor_base[i]);
				}

				if(mpz_cmp_ui(evaluated_poly[j], 1) == 0) {
					++fact_count;
					if(fact_count >= max_fact) {
						go_on = 0;
					}
				}
			}
		}
	}

	mpz_clear(n_root);
	mpz_clear(intermed);

	remove_not_factorized(exponents, evaluated_poly, As, poly_val_num, base_dim);

	// finalize_vector_mpz(&evaluated_poly, poly_val_num);
	// Si noti che questa istruzione apparentemente innocua porta all'uscita di demoni dal naso del programmatore

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
