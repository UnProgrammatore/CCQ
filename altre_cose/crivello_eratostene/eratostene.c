/* Crivello di eratostene parallelo (openmp).
 * Parametri:
 * 	"int * sieve": array di int di lunghezza almeno n/2 opportunamente
 *  	allocato dal chiamante nel quale saranno marchiati con "1" i
 * 		numeri dispari primi. NOTA BENE: "sieve[i] = 1" indica che il
 * 		numero "(i * 2) + 1" è primo (sieve contiene solo numeri dispari)
 *  "long unsigned n": numero fino al quale si desidera calcolare i primi
 *  Compilazione:
 * 		$ gcc -c -lgmp -fopenmp -o eratostene.o eratostene.c
 */

#include <stdlib.h>
#include <omp.h>

void eratosthenes_sieve(unsigned int * sieve, long unsigned n) {
	int n_threads = omp_get_num_threads();
	
	long unsigned chunk = (n/2)/n_threads;

	long unsigned i;
	#pragma omp parallel for schedule(dynamic, chunk)
	for(i = 0; i <= n/2; ++i)
		sieve[i] = 1; 

	for(i = 3; i <= n; i += 2) {
		
		if(i*i > n)
			i = n;  // -> break;

		if(sieve[i/2] == 1) {
			long unsigned j;
			#pragma omp parallel for schedule(dynamic, chunk)
			for(j = i; j <= n/i; j++)
				sieve[(i*j)/2] = 0;

		}
	}
}
