#include "../include/quadratic_sieve.h"

#include <stdio.h>

unsigned long quadratic_sieve(mpz_t N, 
			      unsigned int n, 
			      unsigned int poly_val_num,
			      mpz_t m) {
  double t1, t2;
  
  mpz_t s;
  mpz_init(s);
  mpz_sqrt(s, N); 


  t1 = omp_get_wtime();
  /* Individuazione primi in [2, n] */
  unsigned int * primes = malloc(sizeof(unsigned int) * n);  
  eratosthenes_sieve(primes, n);
  /* Compattiamo i numeri primi in primes */
  unsigned j = 0;
  for(int i = 2; i < n; ++i)
    if(primes[i] == 1) {
      primes[j++] = i;
    }
  unsigned n_all_primes = j;

  unsigned int fattore_semplice = trivial_fact(N, primes, n_all_primes);
  if(fattore_semplice != 0) {
    mpz_set_ui(m, fattore_semplice);
    return OK;
  }

  /* Calcolo base di fattori e soluzioni dell'eq x^2 = N mod p */
  pair * solutions = malloc(sizeof(pair) * n_all_primes);
  unsigned int * factor_base = primes;
    //malloc(sizeof(unsigned int) * n_all_primes); 

  unsigned n_primes = base_fattori(N, s, factor_base, solutions,
				  primes, n_all_primes);
  t2 = omp_get_wtime();
  double t_base = t2 - t1;

  //for(int i = 0; i < n_primes; ++i)
  //  printf("%d\n", factor_base[i]);
 
  printf("dimensione base di fattori: %d\n", n_primes);

 
  /* Vettore degli esponenti in Z */
  t1 = omp_get_wtime();
  unsigned int ** exponents;
  init_matrix(& exponents, poly_val_num, n_primes);
  t2 = omp_get_wtime();
  double t_camp = t2 - t1;
  /* Vettore degli (Ai + s) */
 

  t1 = omp_get_wtime();
  mpz_t * As;
  init_vector_mpz(& As, poly_val_num);

  /* Parte di crivello: troviamo le k+n fattorizzazioni complete */
  unsigned int n_fatt;

 
  n_fatt = sieve(N, factor_base, n_primes, solutions, 
		 exponents, As, poly_val_num, 1000000);
  t2 = omp_get_wtime();
  double t_sieve = t2 - t1;

  printf("numero fattorizzazioni complete trovate: %d\n", n_fatt);

  t1 = omp_get_wtime();
  /* Matrice di esponenti in Z_2 organizzata a blocchi di bit */ 
  word ** M;
  /* Numero di blocchi di bit da utilizzare */
  unsigned long n_blocchi = n_primes / N_BITS + 1;
  /* Inizializzazione egli esponenti mod 2 */
  init_matrix_l(& M, n_fatt, n_blocchi);
  for(int i = 0; i < n_fatt; ++i)
    for(int j = 0; j < n_primes; ++j) {
      unsigned int a = get_matrix(exponents, i, j); 
      set_k_i(M, i, j, a);
    }

  /* Vettore con le info (bit piu' a dx e num bit a 1) su M */
  struct row_stats * wt = malloc(sizeof(struct row_stats) * n_fatt);
  for(int i = 0; i < n_fatt; ++i)
    get_wt_k(M, i, n_primes, & wt[i]);

  /* Eliminazione gaussiana */
  
  gaussian_elimination(exponents, M, As, N, 
		     n_fatt, n_primes, n_blocchi, wt);
  t2 = omp_get_wtime();
  double t_gauss = t2 - t1;

  /* In m ritorno un fattore non banale di N */
  unsigned int n_fact_non_banali = factorization(N, factor_base, 
						 M, exponents, 
						 As, wt, n_fatt, 
						 n_primes, m);

  printf("#time_base time_sieve time_gauss time_totale\n");
  printf("%.6f ", t_base);
  printf("%.6f ", t_sieve);
  printf("%.6f ", t_gauss);
  printf("%.6f ", t_camp);
  //double t_camp = 0;
  printf("%.6f\n", t_base + t_gauss + t_sieve + t_camp);

  if(n_fact_non_banali > 0) {
    /* Pulizia della memoria */
    //finalize_vector(& primes);
    //free(solutions);
    //finalize_matrix(& exponents, poly_val_num);
    finalize_vector_mpz(& As, poly_val_num);
    //finalize_matrix_l(& M, n_fatt);
    free(wt);
    mpz_clear(s);

    return OK;
  }
  else {
    /* Pulizia della memoria */
    //finalize_vector(& primes);
    //free(solutions);
    //finalize_matrix(& exponents, poly_val_num);
    //finalize_vector_mpz(& As, poly_val_num);
    //finalize_matrix_l(& M, n_fatt);
    //free(wt);
    //pz_clear(s);
    
    return SOLO_FATTORIZZAZIONI_BANALI;
  }
}
   
