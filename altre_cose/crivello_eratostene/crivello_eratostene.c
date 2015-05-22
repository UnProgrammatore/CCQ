/* compilation:
 * gcc -fopenmp -lgmp -o crivello_eratostene crivello_eratostene.c
 */

#include <stdlib.h>
#include <omp.h>
#include <gmp.h>

typedef long unsigned lu;

void eratosthenes_sieve(short * sieve, lu n) {
  #pragma opm parallel shared(sieve)
  {

    long unsigned i;
    #pragma omp for schedule(auto)
    for(i = 0; i <= n; i++)
      sieve[i] = 1;
    
    #pragma omp barrier

    mpz_t quad;
    mpz_init(quad);
    mpz_t i_;
    mpz_init(i_);

    #pragma omp for schedule(auto)
    for(i = 2; i <= n; ++i) {
      // uso mpz_t in questo punto per evitare l'overflow
      // quando i*i
      mpz_set_ui(i_, i); 
      mpz_pow_ui(quad, i_ , 2);
      if(mpz_cmp_ui(quad, n) > 0)  // i*i > n
        i = n;
      if(sieve[i] == 1) {
	long unsigned j;
        for(j = i; i*j <= n; j++)
	  sieve[i*j] = 0;
      }
    }
    mpz_clear(quad);
  } // pragma omp parallel
}

int main(int argc, char * argv[]) {
  long unsigned n = atoll(argv[1]);
  short * sieve = malloc(sizeof(short) * n);

#ifdef CRIVELLO_QUAD_VERBOSE
  clock_t c1, c2;
  c1 = clock();
#endif

  eratosthenes_sieve(sieve, n);

#ifdef CRIVELLO_QUAD_VERBOSE 
  c2 = clock();
  printf("#N time\n");
  printf("%lld %f\n", n, (float) (c2 - c1)/CLOCKS_PER_SEC);
#endif

  free(sieve);
}
