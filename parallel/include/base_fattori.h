#include <stdio.h>
#include <gmp.h>

#include "pair.h"
#include "eratostene.h"

short uguali_modulo_mpz_ui(long unsigned f1, mpz_t f2, long unsigned modulo);

pair calcola_soluzioni(mpz_t radice, long unsigned j, long unsigned p);

long unsigned base_fattori(mpz_t numero, mpz_t radice, unsigned int base_fattori[], pair soluzioni[], 
			   unsigned int primi[], long unsigned quanti_primi);
