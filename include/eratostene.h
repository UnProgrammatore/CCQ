/* Crivello di eratostene.
 * Trova i numeri primi in [2, n).
 * Parametri:
 * 	"int * sieve": array di int di lunghezza almeno n opportunamente
 *  	allocato dal chiamante nel quale saranno marchiati con "1" i
 *      primi.
 */

#ifndef ERATOSTENE_H
#define ERATOSTENE_H 1

#include <stdlib.h>
#include <omp.h>

void eratosthenes_sieve(unsigned int * sieve, long unsigned n);

#endif // ERATOSTENE_H
