

#include "factor_base.h"

/*
	Questa funzione deve ritornare una base di fattori per il crivello quadratico.
	Dato N il numero da fattorizzare si considerano i numeri primi, trovati con il 
	crivello di Eratostene minori di exp((1/2)sqrt(log(N)loglog(N))), 
	successivamente vanno esclusi dalla lista quei numeri per i quali il simbolo di 
	Legendre (N|p)!=1.
*/


/*
    Inizializzo all'esterno l'array con il numero giusto secondo l'approssimazione
    di elementi.

    La funzione prende come parametro erat che è un array di dimensione dim_erat inizializzato all'esterno
    mediante il crivello di eratostene, calcola i numeri primi della base p tali per cui
    (N|p)=1 e li mette nel vettore fb di dimensione massima fb_dim.
    La funzione tiene memoria di quanti primi ha inserito in fb, questo valore 
    andrà a sovrascrivere fb_dim.

    num_call è un parametro ausiliario che permette di aumentare la dimensione della base di fattori.

	Scorro tutti gli elementi della base di fattori 
	ritornata dalla funzione eratostene(), 
	quelli che hanno legendre(base.value)!=1 vanno eliminati
	Siamo sicuri che il primo elemento della base ritornata da eratostene sia
	2, quindi non dobbiamo eliminarlo.
*/


/*

*/

void factor_base_erat(mpz_t N, unsigned int* erat, unsigned int dim_erat, unsigned int* fb, unsigned int* fb_dim){

	int res = 0;

	
	for(int i=0; i<dim_erat; ++i){
		int ls = legendre(N,erat[i]);
		if(ls==1){
			fb[res]=erat[i];
			res++;
		}

	}
    
    *fb_dim = res;

}






/*

   
   Precalcolo la dimensione della base di fattori secondo l'approssimazione di Pomerance
	
	mpfr_t pomerance_approx;
	mpfr_init(pomerance_approx);

	double d_exp;
	d_exp = (double)(1/2+num_call)*(log(N)*log(log(N)));
	d_exp = ceil(d_exp);
	unsigned int ui_exp = (unsigned int)d_exp

	mpfr_t e;
	mpfr_init(e);
	mpfr_set_d(e, 2,71828, GMP_RNDN);

    mpfr_pow_ui(pomerance_approx, e, ui_exp, GMP_RNDN );
*/

