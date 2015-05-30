#ifndef GAUSSIAN_ELIMINATION_H
#define GAUSSIAN_ELIMINATION_H 1

/* Tipo di dato per la dimesione del blocco di bit */
#define TYPE unsigned long
/* Numero di bit nel tipo di dato  */
#define N_BITS 64

#include "matrix.h"
#include <stdio.h>

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

/* Funzione che ritorna l'i-mo bit della k-ma
 * riga della matrice M */
unsigned get_k_i(word ** M, unsigned long k, unsigned long i);

/* Funzione che esegue la somma modulo 2 dei vettori
 * v(k) = v(j) + v(i). Utilizzo lo XOR bit a bit che corrisponde
 * alla somma in modulo 2. Eseguo lo XOR tra ogni blocco dei vettori */
unsigned long add_j_to_k(word ** M, unsigned long k, 
			 unsigned long j, unsigned long n_blocks);

/* Funzione che setta la struttura row_stats con le informazioni
 * sulle righe della matrice (ultimo bit a dx e numero bit a 1) */
void get_wt_k(word ** M, unsigned long k, unsigned long n_col, 
	      struct row_stats * wt);

/* Funzione che esegue l'eliminazione gaussiana */
void bit_gaussian_elimination_mod_2(word ** M,
				    unsigned long n_row,
				    unsigned long n_col,
				    unsigned long n_blocks,
				    struct row_stats wt[]);

#endif // GAUSSIAN_ELIMINATION_H
