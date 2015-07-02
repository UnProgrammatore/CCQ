#ifndef SIEVE_GUARD
#define SIEVE_GUARD

#include "pair.h"
#include "vector.h"
#include "matrix.h"
#include "quadratic_sieve.h"

#include <gmp.h>
#include <omp.h>

/* Schema della suddivisione dell'intervallo dei valori di A  
 *
 * begin                        end
 * |______________________________|  size = interval
 *
 * begin_thread   end_thread
 * |______________|                  size = dom_decomp = interval / n_threads
 * 
 *             l  end_block
 * |__||__||__||__|                  size = block_size    
 */

/* Coppia di soluzioni xp e yp in mpz */
struct mpz_pair {
  mpz_t sol1;
  mpz_t sol2;
};
typedef struct mpz_pair mpz_pair;

/* Funzione che esegue la parte di crivello nell'algoritmo:
 * 
 * trova almeno base_dim + max_fact fattorizzazioni completi
 * di elementi Q(A) per A in [begin, begin + interval] e
 * spedisce al main mediante MPI_Send il vettore degli
 * esponenti trovato e (A + s) calcolato. 
 *
 * La funzione ritorna il flag di stop. */
unsigned int
smart_sieve(mpz_t n, // Il numero da fattorizzare
	    unsigned int* factor_base, // La base di fattori
	    unsigned int base_dim, // La dimensione della base di fattori
	    pair* solutions, // Il vettore contenente le soluzioni all'equazione
	    mpz_t begin, // Inizio della sequenza degli A
	    unsigned int interval, // A in [begin, begin + interval]
	    unsigned int block_size, // Porzione di [startfrom, M] da calcolare alla volta
	    unsigned int max_fact); // max_fact + base_dim = massime fattorizzazioni

#endif // SIEVE_GUARD
