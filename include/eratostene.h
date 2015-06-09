/* Crivello di eratostene parallelo (openmp).
 * Parametri:
 * 	"int * sieve": array di int di lunghezza almeno n/2 opportunamente
 *  	allocato dal chiamante nel quale saranno marchiati con "1" i
 * 		numeri dispari primi. NOTA BENE: "sieve[i] = 1" indica che il
 * 		numero "i" è primo (sieve contiene solo numeri dispari)
 *  "long unsigned n": numero fino al quale si desidera calcolare i primi
 *  Compilazione:
 * 		$ gcc -c -lgmp -fopenmp -o eratostene.o eratostene.c
 */

#ifndef ERATOSTENE_H
#define ERATOSTENE_H 1

#include <stdlib.h>
#include <omp.h>

unsigned long eratosthenes_sieve(unsigned int * sieve, long unsigned n);

#endif // ERATOSTENE_H
