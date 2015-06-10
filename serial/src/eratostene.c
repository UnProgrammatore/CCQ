/* Crivello di eratostene parallelo (openmp).
 * Parametri:
 * 	"int * sieve": array di int di lunghezza almeno n/2 opportunamente
 *  	allocato dal chiamante nel quale saranno marchiati con "1" i
 * 	numeri dispari primi. NOTA BENE: "sieve[i] = 1" indica che il
 * 	numero "(i * 2) + 1" è primo (sieve contiene solo numeri dispari)
 *  "long unsigned n": numero fino al quale si desidera calcolare i primi
 *  Compilazione:
 * 		$ gcc -c -lgmp -fopenmp -o eratostene.o eratostene.c
 */

#include "../include/eratostene.h"

void eratosthenes_sieve(unsigned int * sieve, 
				 long unsigned n) {
	unsigned long n_primes = 0;

	for(long unsigned i = 0; i < n; i++)
    		sieve[i] = 1;
	for(long unsigned i = 2; i < n; i++){
		if(i*i > n)
      			break;
    		if(sieve[i] == 1) {
      			for(long unsigned j = i; i*j < n; j++)
				sieve[i*j] = 0;
		}  
	}
}
