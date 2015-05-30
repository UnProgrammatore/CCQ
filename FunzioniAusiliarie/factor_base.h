#include "../include/eratostene.h"
#include <gmp.h>
#include <math.h>
#include <stdlib.h>


/*
	factor_base_erat dà per scontato che qualcuno abbia già calcolato col crivello di eratostene il numero 
	  giusto di numeri primi.
*/
void factor_base_erat(mpz_t N, unsigned int* erat, unsigned int dim_erat, unsigned int* fb, unsigned int& fb_dim);



/*
	factor_base si calcola al suo interno il numero giusto di numeri primi, usando come parametro num_call 
	  che le serve per sapere se aumentare la dimensione della base.
	Al suo interno alloca e dealloca lo spazio necessario per mantenere il vettore di primi derivati 
	  dal calcolo del crivello di eratostene
*/
void factor_base(mpz_t N, unsigned int* fb, unsigned int& fb_dim, unsigned int num_call = 0);