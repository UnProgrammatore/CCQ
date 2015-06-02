#define TYPE unsigned long
#define N_BITS 64
#define MAX_TYPE 18446744073709551615UL

//#include "../include/gaussian_elimination.h"

#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"
#include <omp.h>
#include "sieve.h"
#include "vector.h"
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


/****************/

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

/****************/
/* Funzione che realizza:
 *   a = b * c mod n */
void modular_multiplication(mpz_t a, mpz_t b, mpz_t c, mpz_t n) {
  mpz_mul (a, b, c);
  mpz_mod (a, a, n);
}

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
				mpz_t * Q_A,
				mpz_t N,
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
	//gmp_printf("%Zd * %Zd = ", Q_A[k], Q_A[j]);
	modular_multiplication(Q_A[k], Q_A[k], Q_A[j], N); // Q(Ak) = Q(Ak) * Q(Aj)
	//gmp_printf("%Zd\n", Q_A[k]);
	get_wt_k(M_z2, k, n_col, & wt[k]); // aggiorno wt
      }
    }
    //printf("\n");
    //print_all(M_z2, n_row, n_blocks);
    //printf("\n-----------\n");
  }
}

/* Funzione che ritorna se una riga, nella matrice a blocchi mod 2,
 * è nulla */
int row_is_null(word ** M_z2, unsigned long k, unsigned long n_col, 
		 unsigned long n_blocks) {

  for(unsigned long i = 0; i < n_blocks-1; ++i) {
    if(get_matrix_l(M_z2, k, i) != 0)
      return 0;
  }

  unsigned long n_shift = N_BITS - (n_col % N_BITS);

  /* 01010 ... 110XXX ... X      <- X è il padding 
   *
   * 1111111111111111111111      <- è MAX_TYPE (tutti 1)
   * MAX_TYPE << n_shift =
   * 1111111111111000000000      <- mettendo in & ignoro il padding*/

  print_bits(MAX_TYPE << n_shift);
  printf("\n");
  
  word b = get_matrix_l(M_z2, k, n_blocks-1) & MAX_TYPE << n_shift;
  if(b != 0)
    return 0;
  return 1;  
}

void congruence_relation(mpz_t N, // numero da fattorizzare
			 unsigned int * factor_base,
			 word ** M_z2, // esponenti mod 2
			 unsigned int ** M_z, // esponenti interi
			 mpz_t * Q_a, // vettore Q(Ai)
			 struct row_stats * wt, // zeri sulle righe
			 unsigned long n_row,
			 unsigned long n_primes) {

  mpz_t mpz_temp;
  mpz_init(mpz_temp);

  mpz_t mpz_prime;
  mpz_init(mpz_prime);

  mpz_t X;
  mpz_init(X);

  mpz_t Y;
  mpz_init(Y);

  mpz_t m;
  mpz_init(m);

  mpz_t q;
  mpz_init(q);

  unsigned int exp;


  for(unsigned long i = 0; i < n_row; ++i)
    if(wt[i].n_bit == 0) { // dipendenza trovata
      mpz_set_ui(Y, 1);
      for(int j = 0; j < n_primes; ++j) {
	mpz_set_ui(mpz_prime, factor_base[j]);
	//gmp_printf("prime=%Zd\n", mpz_prime);
	exp = get_matrix(M_z, i, j) / 2;
	//printf("ok\n");
	// temp = (factor_base[j])^(M_z[i][j]) mod N
	mpz_powm_ui(mpz_temp, mpz_prime, exp, N);
	//gmp_printf("temp = %Zd = %Zd^%lu\n", mpz_temp, mpz_prime, exp);
	
	// Y = Y * temp mod N
	modular_multiplication(Y, Y, mpz_temp, N);
	//gmp_printf("Y = %Zd\n", mpz_temp);
      }

      mpz_set(X, Q_a[i]);
      //gmp_printf("(A+s) = %Zd\n", Q_a[i]);
      gmp_printf("mcd(%Zd + %Zd, %Zd) = ", X, Y, N);
      mpz_add(X, X, Y); // X = X + Y
   
      mpz_gcd(m, X, N); // m = mcd(X + Y, N)    
      gmp_printf("%Zd", m);

      mpz_divexact(q, N, m); // q = N / m;

      //gmp_printf("%Zd * %Zd = %Zd, N = ", m, q, N);

      if(mpz_cmp(m, N) < 0 &&  mpz_cmp_ui(m, 1) > 0) { // fatt. non banale
	gmp_printf(", N = %Zd * %Zd\n", m, q);
      }
      else
	printf("\n");
    }

  //mpz_clears(mpz_temp, mpz_prime, X, Y, m, q);
}

/*****************************************************/

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

  print_M(exponents, poly_val_num, n_primes);
  
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
  gaussian_elimination_mod_2(exponents, M, Q_A, N, n_fatt, n_primes, n_blocchi, wt);
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
  */
  //printf("\n");

  congruence_relation(N, factor_base, M, exponents, Q_A, wt, n_fatt, n_primes);

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
