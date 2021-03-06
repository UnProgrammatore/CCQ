#include <omp.h>
#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>

#include "../include/quadratic_sieve.h"

/* usage:
 * ./qs P1 P2 n max_poly_val max_fact interval  
 * 
 * P1 e P2 sono due numeri primi.
 * n è il numero che viene passato al criv. di eratostene
 * max_poly_val è l'M per cui calcoliamo gli A in [0, M]
 * max_fact numero dipendenze linari in più da trovare (dim_base + max_fact)
 * interval è la grandezza di intervalli in cui separare [0, poly_val_num]
 */


int main(int argc, char ** argv) {
  double t1, t2;  

  mpz_t N;
  mpz_init(N);

  mpz_t P1;
  mpz_init(P1);

  mpz_t P2;
  mpz_init(P2);
  
  mpz_set_str(P1, argv[1], 10); 
  mpz_set_str(P2, argv[2], 10); 
  mpz_mul(N, P1, P2);

  gmp_printf("N: %Zd = %Zd * %Zd \n", N, P1, P2);

  unsigned int n = atoi(argv[3]);
  unsigned int poly_val_num = atoi(argv[4]);;
  unsigned int max_fact = atoi(argv[5]);
  unsigned int interval = atoi(argv[6]);
 
  mpz_t N1;
  mpz_t N2;
  mpz_init(N1);
  mpz_init(N2);

  t1 = omp_get_wtime();
  int status = quadratic_sieve(N, n, poly_val_num, max_fact, interval, N1);
  t2 = omp_get_wtime();
  double time = t2 - t1;

  if(status == OK) {
    mpz_divexact(N2, N, N1);
    gmp_printf("#Fattorizzazione trovata: %Zd * %Zd = %Zd\n\n", 
		N1, N2, N);
    
    printf("#N time\n");
    gmp_printf("%Zd %.6f\n", N, time);

    exit(0);
  }
  if(status == SOLO_FATTORIZZAZIONI_BANALI) {
    printf("#Nessuna fattorizzazione non banale trovata\n\n"); 

    //printf("#N time\n");
    //gmp_printf("%Zd %.6f\n", N, time);
    
    exit(1);
  }

}
