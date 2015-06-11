#ifndef LINEAR_ALGEBRA_H
#define LINEAR_ALGEBRA_H 1

/* Tipo di dato per la dimesione del blocco di bit */
#define TYPE unsigned long
/* Numero di bit nel tipo di dato  */
#define N_BITS 64

#include "matrix.h"
#include "vector.h"

#include <stdio.h>
#include <gmp.h>

/* Funzioni per realizzare l'eliminazione gaussiana in modulo 2
 * La matrice di bit degli espondenti modulo 2 sara' organizzata
 * nel seguente modo:
 *
 *       N_BITS        N_BITS            N_BITS
 *    1) [000 ... 001] [000 ... 001] ... [000 ... 001]
 *    2) [000 ... 000] [000 ... 010] ... [000 ... 001]
 *            ...         ...
 *    K) [010 ... 101] [010 ... 101] ... [000 ... 001] */

/* Tipo di dato che determina la dimensione del blocco di bit */
typedef TYPE word;

/* Struttura che contiene le informazioni
 * sulle righe della matrice degli esponenti
 * modulo 2 */
struct row_stats {
  /* posizione ultimo bit piu' a destra */
  long unsigned b_dx;
  /* numero di bit a 1 */
  long unsigned n_bit;
};

/* Funzione che realizza:
 *   a = b * c mod n */
void modular_multiplication(mpz_t a, mpz_t b, mpz_t c, mpz_t n);

/* Funzione che ritorna l'i-mo bit della k-ma
 * riga della matrice M */
unsigned get_k_i(word ** M, unsigned long k, unsigned long i);

/* Funzione che assegna all'i-mo bit della k-ma colonna della
 * matrice a blocchi di bit "M" il valore "value" in mod 2 */
void set_k_i(word ** M, unsigned long k, 
	     unsigned long i, unsigned int value);

/* Funzione che esegue la somma modulo 2 dei vettori
 * v(k) = v(j) + v(k). Utilizzo lo XOR bit a bit che corrisponde
 * alla somma in modulo 2. Eseguo lo XOR tra ogni blocco dei vettori */
void add_vector_z2(word ** M, unsigned long k, 
		   unsigned long j, unsigned long n_blocks);

/* Funzione che esegue la somma in Z dei vettori
 * v(k) = v(j) + v(k) */
void add_vector_z(mpz_t ** M, unsigned long k,
		  unsigned long j, unsigned long n_col);

/* Funzione che setta la struttura row_stats con le informazioni
 * sulle righe della matrice (ultimo bit a dx e numero bit a 1) */
void get_wt_k(word ** M, unsigned long k, unsigned long n_col, 
	      struct row_stats * wt);

/* Funzione che esegue l'eliminazione gaussiana */
void gaussian_elimination(mpz_t ** M_z, // matrice esponenti in Z
			  word ** M_z2, // matrice esponenti in Z2
			  mpz_t * A, // vettore dei (A + s) calcolati
			  mpz_t N, // numero da fattorizzare
			  unsigned long n_row, // #fattorizzaz. complete
			  unsigned long n_col, // numero base di fattori
			  unsigned long n_blocks, // blocchi di bit
			  struct row_stats wt[]); // info sui bit

/* Funzione che presa la matrice ridotta degli esponenti individua 
 * le dipendenze lineari che portano alle congruenze X^2 = Y^2 mod N. 
 * Calcola l'mcd tra (X+Y, N) e restituisce "true" ed una 
 * fattorizzazione non banale ("N1", "N2") di "N" oppure "false" 
 * e N1 = null N2 = null */
unsigned factorization(mpz_t N, // numero da fattorizzare
		       unsigned int * factor_base,
		       word ** M_z2, // esponenti mod 2
		       mpz_t ** M_z, // esponenti interi
		       mpz_t * As, // (Ai + s) moltiplicati tra loro
		       struct row_stats * wt, // zeri sulle righe
		       unsigned long n_row, // #fattorizzaz. complete
		       unsigned long n_primes, // numero base di fattori
		       mpz_t m); // fattore non banale di N

#endif // LINEAR_ALGEBRA_H
