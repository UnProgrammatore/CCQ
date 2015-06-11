#ifndef QUADRATIC_SIEVE_H
#define QUADRATIC_SIEVE_H 1

#include "../include/linear_algebra.h"
#include "../include/vector.h"
#include "../include/matrix.h"
#include "../include/base_fattori.h"
#include "../include/trivial_fact.h"
#include "../include/smart_sieve.h"

#include <omp.h>
#include <gmp.h>
#include <mpi.h>

/* Stati ritornati dall'algoritmo */
enum qs_error_codes{
  OK, // L'algortimo ha trovato una fattorizzazione non banale
  SOLO_FATTORIZZAZIONI_BANALI, // Nessuna fattorizzazione non banale trovata
  NUM_PRIMO, // N su cui si lancia l'algortimo è primo
  IM_A_SLAVE // Ritornato dagli slave
};

unsigned long quadratic_sieve(mpz_t N, // Numero da fattorizzare
			      unsigned int n, // Parametro per il crivello di eratostene
			      unsigned int poly_val_num, // A in [0, poly_val_num]
			      unsigned int max_fact, // Numero di fatt. complete da trovare
			      unsigned int interval, // Numero di intervalli in cui spezzare sieve
			      mpz_t m); // Fattore primo eventualmente trovato

#endif // QUADRATIC_SIEVE_H
