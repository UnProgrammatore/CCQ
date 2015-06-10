#include "../include/quadratic_sieve.h"

#include <stdio.h>

unsigned long quadratic_sieve(mpz_t N, 
			      unsigned int n, 
			      unsigned int poly_val_num,
			      unsigned int max_fact,
			      unsigned int interval,
			      mpz_t m) {
  double t1, t2;
  
  /* Controllo con test di pseudoprimalità di rabin */
  if(mpz_probab_prime_p(N, 25)) {
    return NUM_PRIMO;
  }

  /* Radice intera di N */
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

  unsigned int n_all_primes = j;

  unsigned int simple_factor = trivial_fact(N, primes, n_all_primes);
  if(simple_factor != 0) {
    mpz_set_ui(m, simple_factor);
    return OK;
  }

  /* Calcolo base di fattori e soluzioni dell'eq x^2 = N mod p */
  pair * solutions = malloc(sizeof(pair) * n_all_primes);
  unsigned int * factor_base = primes;

  unsigned n_primes = base_fattori(N, s, factor_base, solutions,
				  primes, n_all_primes);
  t2 = omp_get_wtime();
  double t_base = t2 - t1; 
  printf("dimensione base di fattori: %d\n", n_primes);

  /* Vettore degli esponenti in Z */
  unsigned int ** exponents;
  init_matrix(& exponents, n_primes + max_fact, n_primes);
  t2 = omp_get_wtime();
  double t_camp = t2 - t1;
  /* Vettore degli (Ai + s) */
 
  t1 = omp_get_wtime();
  mpz_t * As;
  init_vector_mpz(& As, n_primes + max_fact);

  /* Parte di crivello: troviamo le k+n fattorizzazioni complete */
  unsigned int n_fatt;
 
  n_fatt = smart_sieve(N, factor_base, n_primes, solutions, 
		 exponents, As, poly_val_num, max_fact, interval);
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

  /* In gauss gli esponenti sommati possono andare in overflow,
     li converto dunque in mpz */
  printf("crasha subito dopo\n");
  mpz_t ** exponents_mpz;
  mpz_t temp;
  mpz_init_set_ui(temp, 2);

  unsigned int a; 
  init_matrix_mpz(& exponents_mpz, n_fatt, n_primes);
  for(unsigned i = 0; i < n_fatt; ++i)
    for(unsigned j = 0; j < n_primes; ++j) {
      a = get_matrix(exponents, i, j);
      mpz_set_ui(temp, a);
      set_matrix_mpz(exponents_mpz, i, j, temp);
    }
  
  printf("gaussian elimination\n");
  /* Eliminazione gaussiana */
  gaussian_elimination(exponents_mpz, M, As, N, 
		       n_fatt, n_primes, n_blocchi, wt);
  t2 = omp_get_wtime();
  double t_gauss = t2 - t1;

  /* In m ritorno un fattore non banale di N */
  unsigned int n_fact_non_banali = factorization(N, factor_base, 
						 M, exponents_mpz, 
						 As, wt, n_fatt, 
						 n_primes, m);
  
  printf("#time_base time_sieve time_gauss time_totale\n");
  printf("%.6f ", t_base);
  printf("%.6f ", t_sieve);
  printf("%.6f ", t_gauss);
  printf("%.6f\n", t_base + t_gauss + t_sieve);
  
  if(n_fact_non_banali > 0) {
    /* Pulizia della memoria */
    return OK;
  }
  else {
    /* Pulizia della memoria */
    return SOLO_FATTORIZZAZIONI_BANALI;
  }
}
   
