#define TYPE unsigned long
#define N_BITS 64

//#include "../include/gaussian_elimination.h"

#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"
#include <omp.h>
#include "sieve.h"
#include <gmp.h>

/* Funzioni per realizzare l'eliminazione gaussiana in modulo 2
 * La matrice di bit degli espondenti modulo 2 sara' organizzata
 * nel seguente modo:
 *
 *       N_BITS        N_BITS            N_BITS    X = eventuale padding
 *    1) [000 ... 001] [000 ... 001] ... [000 ... 0XX]
 *    2) [000 ... 000] [000 ... 010] ... [000 ... 0XX]
 *            ...         ...
 *    K) [010 ... 101] [010 ... 101] ... [000 ... 0XX] */

/* Tipo di dato che determina la dimensione del blocco di bit */
typedef TYPE word;

struct row_stats {
  // bit piu a destra
  long unsigned b_dx;
  // num di bit a 1
  long unsigned n_bit;
};

/* Funzione che ritorna l'i-mo bit della k-ma
 * riga della matrice M */
unsigned get_k_i(word ** M, unsigned long k, 
		     unsigned long i) {
  unsigned long I = i / N_BITS;
  unsigned long n_shift = N_BITS - ((i % N_BITS ) + 1);
  
  return (get_matrix_l(M,k,I) >> n_shift) & 1;
}

void set_k_i(word ** M, unsigned long k, 
	     unsigned long i, unsigned int value) {
  unsigned long I = i / N_BITS;
  unsigned long n_shift = N_BITS - ((i % N_BITS ) + 1);

  word b = get_matrix_l(M, k, I);
  
  //printf("I=%lu, n_s=%lu, ", I, n_shift);
  //print_bits(((unsigned long) value % 2) << n_shift);
  //printf(" - ");

  b = b | (((unsigned long) value % 2UL) << n_shift);

  //printf("b=%lu\n", b);
  //print_bits(b);
  //printf("\n");

  set_matrix_l(M, k, I, b);
}

/* Funzione che esegue la somma modulo 2 dei vettori
 * v(k) = v(j) + v(k). Utilizzo lo XOR bit a bit che corrisponde
 * alla somma in modulo 2. Eseguo lo XOR tra ogni blocco dei vettori */
void add_vector_z2(word ** M, unsigned long k, 
		   unsigned long j, unsigned long n_blocks) {
  for(unsigned long I = 0; I < n_blocks; ++I) {
    word b = get_matrix_l(M, k, I) ^ get_matrix_l(M, j, I);
    set_matrix_l(M, k, I, b); 
  }
}

/* Funzione che esegue la somma in Z dei vettori
 * v(k) = v(j) + v(k) */
void add_vector_z(unsigned int ** M, unsigned long k,
		unsigned long j, unsigned long n_col) {
  for(unsigned long i = 0; i < n_col; ++i) {
    unsigned int sum = get_matrix(M, k, i) + get_matrix(M, j, i);
    set_matrix(M, k, i, sum);
  }
}

/* Funzione che setta la struttura row_stats con le informazioni
 * sulle righe della matrice (ultimo bit a dx e numero bit a 1) */
void get_wt_k(word ** M, unsigned long k, unsigned long n_col, 
	      struct row_stats * wt) {
  // Inizializzo indicando l'ultimo bit nella posizione dopo l'ultima

  //wt->b_dx = n_blocks * N_BITS;

  wt->b_dx = n_col;

  wt->n_bit = 0;

  // Scorro partendo dalla fine fino a trovare il primo 1
  unsigned long i = 0;
  while(get_k_i(M, k, i) == 0 && i < n_col) {
    //printf("%d", get_k_i(M, k, i));
    ++i;
  }

  //printf("\n %lu:%d", i, get_k_i(M, k, i));
  // printf("\n");

  /*
  for(int ii = 0; ii < n_col; ++ii)
    for(int kk=0; kk < 1; kk++)
      printf("%d", get_k_i(M, kk, ii));

  printf("\n");
  */

  // Se ho raggiunto la fine non ci sono bit a 1 ed esco
  if(i >= n_col)
    return;

  wt->b_dx = i;

  for(i = i; i < n_col; ++i)
    if(get_k_i(M, k, i))
      wt->n_bit++;
}

/* Funzione che esegue l'eliminazione gaussiana */
void gaussian_elimination_mod_2(unsigned int ** M_z,
				word ** M_z2,
				unsigned long n_row,
				unsigned long n_col,
				unsigned long n_blocks,
				struct row_stats wt[]) {
 for(unsigned long i = 0; i < n_col; ++i) {
    unsigned long j;
    for(j = 0; j < n_row && wt[j].b_dx != i; ++j)
      ;// avanzo j e basta

    for(unsigned k = j + 1; k < n_row; ++k) {
      if(get_k_i(M_z2, k, i)) { // il bit v(k)(i) deve essere a 1
	add_vector_z2(M_z2, k, j, n_blocks); // v(k) = v(k) + v(j) mod 2
	add_vector_z(M_z, k, j, n_col); // v(k) = v(k) + v(j)
	// moltiplicare i Q(A)
	get_wt_k(M_z2, k, n_col, & wt[k]); // aggiorno wt
      }
    }
  }
}

/*******************************************************/

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
  unsigned long n_fatt = 1000;
  unsigned long n_blocchi = 16;//n_primes / N_BIT

  double t1, t2;

  unsigned int poly_val_num = 5;

  mpz_t N;
  mpz_init(N);

  mpz_set_str(N, "8616460799", 10);

  unsigned int factor_base[15] = {2, 5, 7, 11, 17, 23, 37, 47, 59, 67, 71, 83, 89, 97, 101};

  struct pair solutions[15];
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
 
  unsigned int ** exponents;
  init_matrix(& exponents, 34, n_primes);

  sieve(N, factor_base, n_primes, solutions, exponents, 34);

  print_M(exponents, 34, n_primes);

  t1 = omp_get_wtime();
  init_matrix(& M_z, n_fatt, n_primes);
  init_matrix_l(& M, n_fatt, n_blocchi);
  
  for(int i = 0; i < n_fatt; ++i)
    for(int j = 0; j < n_primes; ++j)
      set_matrix(M_z, i, j, rand() % 10);
 
  for(int i = 0; i < n_fatt; ++i)
    for(int j = 0; j < n_primes; ++j) {
      set_k_i(M, i, j, 0);
    }

  for(int i = 0; i < n_fatt; ++i)
    for(int j = 0; j < n_primes; ++j) {
      unsigned int a = get_matrix(M_z, i, j); 
      set_k_i(M, i, j, a);
    }

  struct row_stats * wt = malloc(sizeof(struct row_stats) * n_fatt);
 
  int n_threads = omp_get_num_threads();
  int chunck = n_fatt/n_threads;

  //#pragma omp parallel for schedule(dynamic, chunck)
  for(int i = 0; i < n_fatt; ++i)
    get_wt_k(M, i, n_primes, & wt[i]);

  t2 = omp_get_wtime();
  double t_set_up = t2 - t1;
  
  t1 = omp_get_wtime();
  gaussian_elimination_mod_2(M_z, M, n_fatt, n_primes, n_blocchi, wt);
  t2 = omp_get_wtime();
  double t_gauss = t2 - t1;

  printf("#time_gauss time_set_up time_totale\n");
  printf("%.6f ", t_gauss);
  printf("%.6f ", t_set_up);
  printf("%.6f\n", t_gauss + t_set_up);
}
