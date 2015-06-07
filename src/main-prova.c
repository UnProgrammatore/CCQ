#include "../include/linear_algebra.h"
#include "../include/vector.h"
#include "../include/sieve.h"
#include "../include/matrix.h"
#include "../include/base_fattori.h"

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <gmp.h>

void print_bits(word a) {
  unsigned int bits[N_BITS];

  for(unsigned int i = 0; i < N_BITS; ++i)
    bits[i] = (a >> i) & 1U;

  for(int i = 63; i >= 0; --i)
    printf("%d", bits[i]);
}

void print_all(unsigned long **M, int righe, int blocchi){
  for(int i = 0; i < righe; ++i) {
    for(int j = 0; j < blocchi; ++j) {
      print_bits(get_matrix_l(M, i, j));
      printf(" ");
    }
    printf("\n");
  }
}

void print_M(unsigned int ** M, int r, int c) {
  for(int i = 0; i < r; ++i) {
   for(int j = 0; j < c; ++j)
     printf("%u, ", get_matrix(M, i, j));
   printf("\n");
  }
}

void print_M_con_i(unsigned int ** M, int r, int c) {
  for(int i = 0; i < r; ++i) {
    printf("%d:  ", i);
   for(int j = 0; j < c; ++j)
     printf("%u, ", get_matrix(M, i, j));
   printf("\n");
  }
}

void print_M_2(unsigned int ** M, int r, int c) {
  for(int i = 0; i < r; ++i) {
   for(int j = 0; j < c; ++j)
     printf("%u", get_matrix(M, i, j) % 2);
   printf("\n");
  }
}


int main() {
  double t1, t2;  

  mpz_t N;
  mpz_init(N);
  // mpz_set_str(N, "8616460799", 10); 
  // mpz_set_str(N, "276833100228154273", 10); 19 0.79

  mpz_t P1;
  mpz_init(P1);

  mpz_t P2;
  mpz_init(P2);

  //mpz_set_str(P1, "1888888888", 10); 
  //mpz_set_str(P2, "5245978939", 10); 

  mpz_set_str(P1, "1888888888", 10); 
  mpz_set_str(P2, "5245978939", 10); 	
  mpz_mul(N, P1, P2);

  gmp_printf("N: %Zd = %Zd * %Zd \n", N, P1, P2);

  //mpz_set_str(N, "439389214701485110197221", 10); 
  mpz_set_str(N, "304285060515912", 10); 

  mpz_t s;
  mpz_init(s);
  mpz_sqrt(s, N);

  unsigned int n = 1000;
  unsigned int * numbers = malloc(sizeof(unsigned int) * n);
  unsigned n_all_primes = eratosthenes_sieve(numbers, n);

  unsigned int * primi = malloc(sizeof(unsigned int) * n_all_primes);
  unsigned int * factor_base = primi;//malloc(sizeof(unsigned int) * n_all_primes);

  unsigned j = 0;
  for(int i = 2; i < n; ++i)
    if(numbers[i] == 1)
      primi[j++] = i;

  pair * solutions = malloc(sizeof(pair) * n_all_primes);

  t1 = omp_get_wtime();
  unsigned n_primes = base_fattori(N, s, factor_base, solutions,
				   primi, n_all_primes);
  t2 = omp_get_wtime();
  double t_base = t2 - t1;

  printf("dimensione base di fattori: %d\n", n_primes);

  /*
  for(int i = 0; i < n_primes; ++i)
    printf("%d\n", factor_base[i]);
  */

  unsigned int poly_val_num = 20000;//12800;

  unsigned int ** exponents;
  init_matrix(& exponents, poly_val_num, n_primes);
  for(int i = 0; i < poly_val_num; ++i)
    for(int j = 0; j < n_primes; ++j)
    set_matrix(exponents, i, j, 0);
  
  mpz_t * As;
  init_vector_mpz(& As, poly_val_num);

  unsigned int n_fatt;
  t1 = omp_get_wtime();
  n_fatt = sieve(N, factor_base, n_primes, solutions, 
		 exponents, As, poly_val_num, 40);
  t2 = omp_get_wtime();
  double t_sieve = t2 - t1;

  printf("numero fattorizzazioni complete trovate: %d\n", n_fatt);

  word ** M;
  unsigned long n_blocchi = n_primes / N_BITS + 1;
  init_matrix_l(& M, n_fatt, n_blocchi);
  /*
  for(int i = 0; i < n_fatt; ++i)
    for(int j = 0; j < n_primes; ++j) {
      set_k_i(M, i, j, 0);
    }
  */
  for(int i = 0; i < n_fatt; ++i)
    for(int j = 0; j < n_primes; ++j) {
      unsigned int a = get_matrix(exponents, i, j); 
      set_k_i(M, i, j, a);
    }

  struct row_stats * wt = malloc(sizeof(struct row_stats) * n_fatt);
  for(int i = 0; i < n_fatt; ++i)
    get_wt_k(M, i, n_primes, & wt[i]);

  t1 = omp_get_wtime();
  gaussian_elimination(exponents, M, As, N, 
		       n_fatt, n_primes, n_blocchi, wt);
  t2 = omp_get_wtime();
  double t_gauss = t2 - t1;

  mpz_t temp;
  mpz_init(temp);
  mpz_t N1;
  mpz_t N2;
  mpz_init(N1);
  mpz_init(N2);
  unsigned int n_fact_non_banali = factorization(N, factor_base, 
						 M, exponents, 
						 As, wt, n_fatt, 
						 n_primes, N1);
  if(n_fact_non_banali > 0) {
    mpz_divexact(N2, N, N1);
    gmp_printf("Fattorizzazione trovata: %Zd * %Zd = %Zd\n\n", 
	       N1, N2, N);
    printf("Numero di fattorizzazioni non banali: %d\n", 
	   n_fact_non_banali);
  }
  else
    printf("Nessuna fattorizzazione non banale trovata\n\n");
  
  printf("#time_base time_sieve time_gauss time_totale\n");
  printf("%.6f ", t_base);
  printf("%.6f ", t_sieve);
  printf("%.6f ", t_gauss);
  printf("%.6f\n", t_base + t_gauss + t_sieve);
}
