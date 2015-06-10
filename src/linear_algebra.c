#include "../include/linear_algebra.h"

#include <omp.h>

void modular_multiplication(mpz_t a, mpz_t b, mpz_t c, mpz_t n) {
  mpz_mul (a, b, c);
  mpz_mod (a, a, n);
}

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

  b = b | (((unsigned long) value % 2UL) << n_shift);

  set_matrix_l(M, k, I, b);
}

void add_vector_z2(word ** M, unsigned long k, 
		   unsigned long j, unsigned long n_blocks) {
  for(unsigned long I = 0; I < n_blocks; ++I) {
    word b = get_matrix_l(M, k, I) ^ get_matrix_l(M, j, I);
    set_matrix_l(M, k, I, b); 
  }
}

void add_vector_z(mpz_t ** M, unsigned long k,
		unsigned long j, unsigned long n_col) {
  mpz_t sum;
  mpz_init(sum);

  mpz_t x;
  mpz_init(x);

  mpz_t y;
  mpz_init(y);
  
  for(unsigned long i = 0; i < n_col; ++i) {
    get_matrix_mpz(x, M, k, i);
    get_matrix_mpz(y, M, j, i);

    mpz_add(sum, x, y); // M[k][i] = M[k][i] + M[j][i]       
    
    set_matrix_mpz(M, k, i, sum);
  }
  
  mpz_clears(sum, x, y, NULL);
}

void get_wt_k(word ** M, unsigned long k, unsigned long n_col, 
	      struct row_stats * wt) {
  // Inizializzo indicando l'ultimo bit nella posizione dopo l'ultima
  wt->b_dx = n_col;
  // Numero bit a 1 = 0
  wt->n_bit = 0;

  // Scorro partendo dalla fine fino a trovare il primo 1
  unsigned long i = 0;
  while(get_k_i(M, k, i) == 0 && i < n_col)
    ++i;

  // Se ho raggiunto la fine non ci sono bit a 1 ed esco
  if(i >= n_col)
    return;

  wt->b_dx = i;

  for(i = i; i < n_col; ++i)
    if(get_k_i(M, k, i))
      wt->n_bit++;
}

void gaussian_elimination(mpz_t ** M_z,
			  word ** M_z2,
			  mpz_t * As,
			  mpz_t N,
			  unsigned long n_row,
			  unsigned long n_col,
			  unsigned long n_blocks,
			  struct row_stats wt[]) {
  double t1, t2;
  double t_Z2 = 0;
  double t_Z = 0;
  double t_As = 0;
  double t_get_wt = 0;
  
  for(unsigned long i = 0; i < n_col; ++i) {
    unsigned long j;
    for(j = 0; j < n_row && wt[j].b_dx != i; ++j)
      ; // avanzo j e basta

    for(unsigned k = j + 1; k < n_row; ++k) {
      
      if(get_k_i(M_z2, k, i)) { // il bit v(k)(i) deve essere a 1
	add_vector_z2(M_z2, k, j, n_blocks); // v(k) = v(k) + v(j) mod 2
     	add_vector_z(M_z, k, j, n_col); // v(k) = v(k) + v(j)
	// (A_k + s) = (A_k + s) * (A_j + s)
	modular_multiplication(As[k], As[k], As[j], N); 
	get_wt_k(M_z2, k, n_col, & wt[k]); // aggiorno wt
      }
    }
  }	
}

unsigned factorization(mpz_t N, // numero da fattorizzare
		       unsigned int * factor_base,
		       word ** M_z2, // esponenti mod 2
		       mpz_t ** M_z, // esponenti interi
		       mpz_t * As, // (Ai + s) moltiplicati tra loro
		       struct row_stats * wt, // zeri sulle righe
		       unsigned long n_row, // #fattorizzaz. complete
		       unsigned long n_primes, // numero base di fattori
		       mpz_t m) { // fattore non banale di N

  mpz_t mpz_temp;
  mpz_init(mpz_temp);

  mpz_t mpz_prime;
  mpz_init(mpz_prime);

  mpz_t X;
  mpz_init(X);

  mpz_t Y;
  mpz_init(Y);

  mpz_t q;
  mpz_init(q);

  mpz_t exp;
  mpz_init(exp);

  unsigned int n_dip = 0;

  for(unsigned long i = 0; i < n_row; ++i)
    if(wt[i].n_bit == 0) { // dipendenza trovata

      ++n_dip;

      mpz_set_ui(Y, 1);
      for(int j = 0; j < n_primes; ++j) {
	mpz_set_ui(mpz_prime, factor_base[j]);
	get_matrix_mpz(exp, M_z, i, j);
	mpz_divexact_ui(exp, exp, 2); // exp = exp / 2
	// temp = (factor_base[j])^(M_z[i][j]) mod N
	mpz_powm(mpz_temp, mpz_prime, exp, N);
	// Y = Y * temp mod N
	modular_multiplication(Y, Y, mpz_temp, N);
      }

      mpz_set(X, As[i]);
      mpz_add(X, X, Y); // X = X + Y   
      mpz_gcd(m, X, N); // m = mcd(X + Y, N)  
     
      mpz_divexact(q, N, m); // q = N / m;
      if(mpz_cmp(m, N) < 0 && mpz_cmp_ui(m, 1) > 0) {
	return 1;
      }
    }

  printf("n_dip = %d\n", n_dip);

  mpz_clears(mpz_temp, exp, q, mpz_prime, X, Y, NULL);

  return n_fatt_non_banali;
}
