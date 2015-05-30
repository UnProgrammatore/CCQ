#include "../include/gaussian_elimination.h"
#include <stdlib.h>

/* Funzioni per realizzare l'eliminazione gaussiana in modulo 2
 * La matrice di bit degli espondenti modulo 2 sara' organizzata
 * nel seguente modo:
 *
 *       N_BITS        N_BITS            N_BITS
 *    1) [000 ... 001] [000 ... 001] ... [000 ... 001]
 *    2) [000 ... 000] [000 ... 010] ... [000 ... 001]
 *            ...         ...
 *    K) [010 ... 101] [010 ... 101] ... [000 ... 001] */

/* Funzione che ritorna l'i-mo bit della k-ma
 * riga della matrice M */
unsigned get_k_i(word ** M, unsigned long k, 
		     unsigned long i) {
  unsigned long I = i / N_BITS;
  unsigned long n_shift = N_BITS - ((i % N_BITS ) + 1);
  
  return (get_matrix_l(M,k,I) >> n_shift) & 1;
}

/* Funzione che esegue la somma modulo 2 dei vettori
 * v(k) = v(j) + v(i). Utilizzo lo XOR bit a bit che corrisponde
 * alla somma in modulo 2. Eseguo lo XOR tra ogni blocco dei vettori */
unsigned long add_j_to_k(word ** M, unsigned long k, 
			 unsigned long j, unsigned long n_blocks) {
  for(unsigned long I = 0; I < n_blocks; ++I) {
    word b = get_matrix_l(M, k, I) ^ get_matrix_l(M, j, I);
    set_matrix_l(M, k, I, b); 
  }
}

/* Funzione che setta la struttura row_stats con le informazioni
 * sulle righe della matrice (ultimo bit a dx e numero bit a 1) */
void get_wt_k(word ** M, unsigned long k, unsigned long n_col, 
	      struct row_stats * wt) {
  // Inizializzo indicando l'ultimo bit nella posizione dopo l'ultima
  wt->b_dx = n_col;

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

/* Funzione che esegue l'eliminazione gaussiana */
void bit_gaussian_elimination_mod_2(unsigned long ** M,
				    unsigned long n_row,
				    unsigned long n_col,
				    unsigned long n_blocks,
				    struct row_stats wt[]) {
 for(unsigned long i = 0; i < n_col; ++i) {
    unsigned long j;
    for(j = 0; j < n_row && wt[j].b_dx != i; ++j)
      ;// avanzo j e basta

    for(unsigned k = j + 1; k < n_row; ++k) {
      if(get_k_i(M, k, i)) { // il bit v(k)(i) deve essere a 1
	add_j_to_k(M, k, j, n_blocks); // v(k) = v(k) + v(j)
	// sommare le righe della matrice degli esponenti in Z
	// moltiplicare i Q(A)
	get_wt_k(M, k, n_col, & wt[k]); // aggiorno wt
      }
    }
  }
}
