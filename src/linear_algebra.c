#include "../include/linear_algebra.h"

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

void add_vector_z(unsigned int ** M, unsigned long k,
		unsigned long j, unsigned long n_col) {
  for(unsigned long i = 0; i < n_col; ++i) {
    unsigned int sum = get_matrix(M, k, i) + get_matrix(M, j, i);
    set_matrix(M, k, i, sum);
  }
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

void gaussian_elimination(unsigned int ** M_z,
				word ** M_z2,
				mpz_t * As,
				mpz_t N,
				unsigned long n_row,
				unsigned long n_col,
				unsigned long n_blocks,
				struct row_stats wt[]) {

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
		       unsigned int ** M_z, // esponenti interi
		       mpz_t * As, // (Ai + s) moltiplicati tra loro
		       struct row_stats * wt, // zeri sulle righe
		       unsigned long n_row, // #fattorizzaz. complete
		       unsigned long n_primes, // numero base di fattori
		       mpz_t m, // fattore non banale di N
		       mpz_t q) { // altro fattore non banale di N

  mpz_t mpz_temp;
  mpz_init(mpz_temp);

  mpz_t mpz_prime;
  mpz_init(mpz_prime);

  mpz_t X;
  mpz_init(X);

  mpz_t Y;
  mpz_init(Y);

  unsigned int exp;

  unsigned int n_dip = 0;
  unsigned int n_fatt_non_banali = 0;

  for(unsigned long i = 0; i < n_row; ++i)
    if(wt[i].n_bit == 0) { // dipendenza trovata

      ++n_dip;

      mpz_set_ui(Y, 1);
      for(int j = 0; j < n_primes; ++j) {
	mpz_set_ui(mpz_prime, factor_base[j]);
	exp = get_matrix(M_z, i, j) / 2;
	// temp = (factor_base[j])^(M_z[i][j]) mod N
	mpz_powm_ui(mpz_temp, mpz_prime, exp, N);
	// Y = Y * temp mod N
	modular_multiplication(Y, Y, mpz_temp, N);
      }

      mpz_set(X, As[i]);
      mpz_add(X, X, Y); // X = X + Y
   
      mpz_gcd(m, X, N); // m = mcd(X + Y, N)    

      mpz_divexact(q, N, m); // q = N / m;

      gmp_printf("%Zd * %Zd\n", m, q);

      if(mpz_cmp(m, N) < 0 &&  mpz_cmp_ui(m, 1) > 0) { // fatt. non banale
	++n_fatt_non_banali;
	return 1;
      }
    }
  return 0;

  mpz_clears(mpz_temp, mpz_prime, X, Y, NULL);
}
