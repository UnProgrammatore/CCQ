#include "../include/linear_algebra.h"
#include "../include/vector.h"
#include "../include/sieve.h"
#include "../include/matrix.h"

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
  word ** M;
  unsigned int ** M_z;

  unsigned long n_primes = 15;
  unsigned long n_blocchi = 1;//n_primes / N_BIT

  double t1, t2;

  unsigned int poly_val_num = 12800;

  mpz_t N;
  mpz_init(N);
 
  mpz_set_str(N, "8616460799", 10);

  unsigned int factor_base[15] = {2, 5, 7, 11, 17, 23, 37, 47, 59, 67, 71, 83, 89, 97, 101};

  pair solutions[15];
  unsigned c = 0;
  solutions[c].sol1 = 1;
  solutions[c++].sol2 = 1;
  solutions[c].sol1 = 3;
  solutions[c++].sol2 = 4;
  solutions[c].sol1 = 6;
  solutions[c++].sol2 = 0;
  solutions[c].sol1 = 6;
  solutions[c++].sol2 = 4;
  solutions[c].sol1 = 15;
  solutions[c++].sol2 = 11;
  solutions[c].sol1 = 7;
  solutions[c++].sol2 = 1;
  solutions[c].sol1 = 21;
  solutions[c++].sol2 = 34;
  solutions[c].sol1 = 15;
  solutions[c++].sol2 = 34;
  solutions[c].sol1 = 9;
  solutions[c++].sol2 = 16;
  solutions[c].sol1 = 51;
  solutions[c++].sol2 = 25;
  solutions[c].sol1 = 69;
  solutions[c++].sol2 = 19;
  solutions[c].sol1 = 68;
  solutions[c++].sol2 = 38;
  solutions[c].sol1 = 8;
  solutions[c++].sol2 = 87;
  solutions[c].sol1 = 52;
  solutions[c++].sol2 = 55;
  solutions[c].sol1 = 34;
  solutions[c++].sol2 = 57;
  /*
  for(int k = 0; k < n_primes; ++k)
    printf("%d: xp=%d, yp=%d\n", factor_base[k], solutions[k].sol1, solutions[k].sol2);
  printf("\n");
  */

  t1 = omp_get_wtime();
 
  unsigned int ** exponents;
  init_matrix(& exponents, poly_val_num, n_primes);
  for(int i = 0; i < poly_val_num; ++i)
    for(int j = 0; j < n_primes; ++j)
      set_matrix(exponents, i, j, 0);

  //print_M(exponents, poly_val_num, n_primes);
  
  mpz_t * Q_A;
  init_vector_mpz(& Q_A, poly_val_num);

  unsigned int n_fatt;
  n_fatt = sieve(N, factor_base, n_primes, solutions, exponents, Q_A, poly_val_num);
  //printf("\n");
  //printf("n_fatt:%d\n\n", n_fatt);

  //print_M_con_i(exponents, poly_val_num, n_primes);

  //init_matrix(& M_z, n_fatt, n_primes);
  init_matrix_l(& M, n_fatt, n_blocchi);

  //unsigned int f_c[] = {34, 453, 1134, 3143, 3388, 4514, 4808, 5251, 6033, 6263, 6683, 7508, 8494, 9086, 10233, 12379, 12799};

  //for(int i = 0; i < n_fatt; ++i)
  //  for(int j = 0; j < n_primes; ++j)
  //    set_matrix(M_z, i, j, get_matrix(exponents, f_c[i], j));

  //print_M(exponents, n_fatt, n_primes);
  
  /*
  for(int i = 0; i < n_fatt; ++i)
    for(int j = 0; j < n_primes; ++j)
      set_matrix(M_z, i, j, rand() % 10);
  */

  for(int i = 0; i < n_fatt; ++i)
    for(int j = 0; j < n_primes; ++j) {
      set_k_i(M, i, j, 0);
    }
 
  for(int i = 0; i < n_fatt; ++i)
    for(int j = 0; j < n_primes; ++j) {
      unsigned int a = get_matrix(exponents, i, j); 
      set_k_i(M, i, j, a);
    }

  //printf("\n");

  //print_all(M, n_fatt, n_blocchi);

  struct row_stats * wt = malloc(sizeof(struct row_stats) * n_fatt);
 
  int n_threads = omp_get_num_threads();
  int chunck = n_fatt/n_threads;

  //#pragma omp parallel for schedule(dynamic, chunck)
  for(int i = 0; i < n_fatt; ++i)
    get_wt_k(M, i, n_primes, & wt[i]);

  t2 = omp_get_wtime();
  double t_set_up = t2 - t1;
  
  t1 = omp_get_wtime();
  gaussian_elimination(exponents, M, Q_A, N, n_fatt, n_primes, n_blocchi, wt);
  t2 = omp_get_wtime();
  double t_gauss = t2 - t1;

  //printf("\n\ngauss:\n");

  //print_M(exponents, n_fatt, n_primes);

  //printf("\n");

  //print_all(M, n_fatt, n_blocchi);

  //printf("\n");

  mpz_t temp;
  mpz_init(temp);
  /*
  for(int i = 0; i < n_fatt; ++i) {
    //mpz_sqrt(temp, Q_A[i]);
    //gmp_printf ("%Zd\n", temp);
    gmp_printf ("%Zd\n", Q_A[i]);
  }
  printf("\n");
  */

  mpz_t N1;
  mpz_t N2;

  mpz_init(N1);
  mpz_init(N2);

  if(factorization(N, factor_base, M, exponents, Q_A, wt, n_fatt, n_primes, N1, N2))
    gmp_printf ("Fattorizzazione trovata: %Zd * %Zd = %Zd\n\n", N1, N2, N);
  else
    printf("Nessuna fattorizzazione non banale trovata\n\n");

  /*
  if(row_is_null(M, 16, n_primes, n_blocchi))
    printf("test1 ok\n");
  else
    printf("test1 errore\n");

  set_matrix_l(M, 16, 0, 1);

  print_bits(get_matrix_l(M, 16, 0));
  printf("\n");

  if(row_is_null(M, 16, n_primes, n_blocchi))
    printf("test2 ok\n");
  else
    printf("test2 errore\n");
  */

  printf("#time_gauss time_set_up time_totale\n");
  printf("%.6f ", t_gauss);
  printf("%.6f ", t_set_up);
  printf("%.6f\n", t_gauss + t_set_up);
}
