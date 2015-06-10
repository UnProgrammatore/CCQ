/* Crivello di eratostene.
 * Trova i numeri primi in [2, n).
 * Parametri:
 * 	"int * sieve": array di int di lunghezza almeno n opportunamente
 *  	allocato dal chiamante nel quale saranno marchiati con "1" i primi.
 * 	"long unsigned n": numero fino al quale si desidera calcolare i primi
 * Ritororna i numeri dispari primi fino ad n. NOTA BENE: "sieve[i] = 1" indica che il
 * numero "i" è primo 
 * Compilazione:
 * 		$ gcc -c -lgmp -fopenmp -o eratostene.o eratostene.c
>>>>>>> 60b860d79450927c9a8fb282afcc0decb489e370
 */

#ifndef ERATOSTENE_H
#define ERATOSTENE_H 1

#include <stdlib.h>
#include <omp.h>

void eratosthenes_sieve(unsigned int * sieve, long unsigned n);

#endif // ERATOSTENE_H
