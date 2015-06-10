#include "../include/smart_sieve.h"
#include "../include/linear_algebra.h"
#include "../include/matrix.h"

/*#include <stdio.h>

void print_M(unsigned int ** M, int r, int c) {
  for(int i = 0; i < r; ++i) {
   for(int j = 0; j < c; ++j)
     printf("%u, ", get_matrix(M, i, j));
   printf("\n");
  }
}
*/

unsigned int smart_sieve(
	mpz_t n,
	unsigned int* factor_base,
	unsigned int base_dim,
	pair* solutions,
	unsigned int** exponents,
	mpz_t* As,
	unsigned int poly_val_num,
	unsigned int max_fact,
	unsigned int intervals
	) {

	unsigned int i; // Indice generico

	unsigned char go_on = 1;

	mpz_t piu_intermed;
	mpz_init(piu_intermed);

	mpz_t ancora_piu_intermed;
	mpz_init(ancora_piu_intermed);

	unsigned int h; // Usato per copiare gli base_dim esponenti

	unsigned int** expo2; // Matrice temporanea degli esponenti
	init_matrix(&expo2, intervals, base_dim);
	word** is_used_expo2; // Vettore che segna quali esponenti sono stati inizializzati
	init_matrix_l(&is_used_expo2, 1, (intervals / N_BITS) + 1);

	mpz_t n_root;
	mpz_t intermed; // Valore appoggio
	mpz_init(n_root);
	mpz_init(intermed);
	mpz_sqrt(n_root, n);

	unsigned int fact_count = 0; // Numero fattorizzazioni trovate
	unsigned int j, k, l; // Indici generici
	mpz_t* evaluated_poly; // Vettore temporaneo dei Q(A) valutati
	mpz_t* As2; // A + s temporanei
	init_vector_mpz(&As2, intervals);

	init_vector_mpz(&evaluated_poly, intervals);

	max_fact += base_dim; // k + n

	for(l = 0; l <= poly_val_num - intervals && go_on; l += intervals) {
		//printf("Giro con l = %d\n", l);
		for(i = 0; i < ((intervals / N_BITS) + 1); ++i) {
			set_matrix_l(is_used_expo2, 0, i, 0);
		}

		for(i = 0; i < intervals; ++i) {
			mpz_add_ui(intermed, n_root, i + l);
			mpz_mul(intermed, intermed, intermed);
			mpz_sub(evaluated_poly[i], intermed, n);
			
			mpz_add_ui(As2[i], n_root, i + l);
		}

		for(i = 0; i < base_dim && go_on; ++i) {
			for(j = solutions[i].sol1; j < intervals + l && go_on; j += factor_base[i]) {

				//printf("Sol1 Sto per provare a fattorizzare %d col primo %d\n", j, factor_base[i]);
				while(mpz_divisible_ui_p(evaluated_poly[j - l], factor_base[i])) {
					
					//gmp_printf("%Zd è divisibile per %d\n", evaluated_poly[j - l], factor_base[i]);
					// Se non sono mai stati usati gli esponenti
					if(get_k_i(is_used_expo2, 0, j - l) == 0) {
						//printf("Quasi al for\n");
						for(k = 0; k < base_dim; ++k)
							set_matrix(expo2, j - l, k, 0);
						//printf("Fatto il for\n");
						set_k_i(is_used_expo2, 0, j - l, 1);
						//printf("Ho fatto i paciughi per il coso\n");
					}

					set_matrix(expo2, j - l, i, get_matrix(expo2, j - l, i) + 1); // ++exponents[j][i];
					mpz_divexact_ui(evaluated_poly[j - l], evaluated_poly[j - l], factor_base[i]);
				}
				if(mpz_cmp_ui(evaluated_poly[j - l], 1) == 0) {
					//printf("Sto per aggiungere il %desimo numero a exponents e l'altro\n", fact_count);
					for(h = 0; h < base_dim; ++h) {
						set_matrix(exponents, fact_count, h, get_matrix(expo2, j - l, h));
					}
					mpz_set(As[fact_count], As2[j - l]);
					unsigned int q;
					mpz_set_ui(intermed, 1);
					for(q = 0; q < base_dim; ++q) {
						mpz_set_ui(ancora_piu_intermed, factor_base[q]);
						mpz_pow_ui(piu_intermed, ancora_piu_intermed, get_matrix(expo2, j - l, q));
						//gmp_printf("AHAHAHAHAH: %Zd\n", piu_intermed);
						mpz_mul(intermed, intermed, piu_intermed);
					}
					//gmp_printf("Ho fattorizzato il numero %Zd in ", intermed);
					for(q = 0; q < base_dim; ++q) {
						//printf("%d^%d ", factor_base[q], get_matrix(expo2, j - l, q));
					}
					//printf("\n");
					++fact_count;
					if(fact_count >= max_fact)
						go_on = 0;
					//printf("Fatto\n");
				}
			}
			solutions[i].sol1 = j; // Al prossimo giro ricominciamo da dove abbiamo finito

			for(j = solutions[i].sol2; j < intervals + l && factor_base[i] != 2 && go_on; j += factor_base[i]) {

				//printf("Sol2 Sto per provare a fattorizzare %d col primo %d\n", j, factor_base[i]);
				while(mpz_divisible_ui_p(evaluated_poly[j - l], factor_base[i])) {
					
					//gmp_printf("%Zd è divisibile per %d\n", evaluated_poly[j - l], factor_base[i]);
					// Se non sono mai stati usati gli esponenti
					if(get_k_i(is_used_expo2, 0, j - l) == 0) {
						//printf("Quasi al for\n");
						for(k = 0; k < base_dim; ++k)
							set_matrix(expo2, j - l, k, 0);
						//printf("Fatto il for\n");
						set_k_i(is_used_expo2, 0, j - l, 1);
						//printf("Ho fatto i paciughi per il coso\n");
					}

					set_matrix(expo2, j - l, i, get_matrix(expo2, j - l, i) + 1); // ++exponents[j][i];
					mpz_divexact_ui(evaluated_poly[j - l], evaluated_poly[j - l], factor_base[i]);
				}
				if(mpz_cmp_ui(evaluated_poly[j - l], 1) == 0) {
					//printf("Sto per aggiungere il %desimo numero a exponents e l'altro\n", fact_count);
					for(h = 0; h < base_dim; ++h) {
						set_matrix(exponents, fact_count, h, get_matrix(expo2, j - l, h));
					}
					mpz_set(As[fact_count], As2[j - l]);
					unsigned int q;
					mpz_set_ui(intermed, 1);
					for(q = 0; q < base_dim; ++q) {
						mpz_set_ui(ancora_piu_intermed, factor_base[q]);
						mpz_pow_ui(piu_intermed, ancora_piu_intermed, get_matrix(expo2, j - l, q));
						//gmp_printf("AHAHAHAHAH: %Zd\n", piu_intermed);
						mpz_mul(intermed, intermed, piu_intermed);
					}
					//gmp_printf("Ho fattorizzato il numero %Zd in ", intermed);
					for(q = 0; q < base_dim; ++q) {
						//printf("%d^%d ", factor_base[q], get_matrix(expo2, j - l, q));
					}
					//printf("\n");
					++fact_count;
					if(fact_count >= max_fact)
						go_on = 0;
					//printf("Fatto\n");
				}
			}
			solutions[i].sol2 = j; // Al prossimo giro ricominciamo da dove abbiamo finito
		}
	}
	//print_M(exponents, fact_count, base_dim);
	return fact_count;
}