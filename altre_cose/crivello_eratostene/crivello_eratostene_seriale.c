/* compilation:
 * gcc crivello_eratostene.c -fopenmp -lgmp -o crivello_eratostene
 */

#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <gmp.h>
#include <time.h>


void eratosthenes_sieve(int * sieve, long unsigned n) {

	long unsigned i;
	for(i = 0; i <= n/2; ++i)
		sieve[i] = 1; 

	//#pragma omp parallel for schedule(dynamic)
	for(i = 3; i <= n; i += 2) {
		
		if(i*i > n)
			i = n;

		if(sieve[i/2] == 1) {
			long unsigned j;
			for(j = i; j <= n/i; j++)
				sieve[(i*j)/2] = 0;

		}
	}
}

int main(int argc, char * argv[]) {
	long unsigned n = atoll(argv[1]);
	int * sieve = (int *) malloc(sizeof(int) * (n/2));

#ifdef TIME
	clock_t c1, c2;
	c1 = clock();
#endif

	eratosthenes_sieve(sieve, n);

#ifdef TIME 
	c2 = clock();
	printf("#N time\n");
	printf("%ld %f\n", n, (float) (c2 - c1)/CLOCKS_PER_SEC);
#endif
#ifdef VERBOSE
	printf("#");
	long unsigned i;
	for(i = 2; i < n+1; ++i)
		if(sieve[i] == 1)
			printf("%lu ", i);
	printf("\n");
#endif
	//free(sieve);

	return 0;
}
