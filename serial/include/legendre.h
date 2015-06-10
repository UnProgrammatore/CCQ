#ifndef LEGENDRE_GUARD
#define LEGENDRE_GUARD

#include <gmp.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>



typedef
struct pair_{
	unsigned int sol1;
	unsigned int sol2;
} pair;


void inizializza_pair(pair p);


/*
	semplice calcolo del simbolo di Legendre:

	    	{	1 se p|/a e a è un residuo quadratico modulo p
	(a|p) = |	0 se p|a
	    	{   -1 se p|/a e a non è un residuo quadratico modulo p
*/
int legendre(mpz_t n, unsigned int p );


/*
	Oltre al calcolo ritorna le soluzioni dell'equazione x^2 = a mod p in 
	pair.

*/
int legendre_sol(mpz_t n, unsigned int p, pair sol );

#endif // LEGENDRE_GUARD