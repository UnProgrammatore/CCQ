#include <gmp.h>
#include<stdlib.h>


/*
	Questa funzione deve ritornare una base di fattori per il crivello quadratico.
	Dato N il numero da fattorizzare si considerano i numeri primi, trovati con il 
	crivello di Eratostene minori di exp((1/2)sqrt(log(N)loglog(N))), 
	successivamente vanno esclusi dalla lista quei numeri per i quali il simbolo di 
	Legendre (N|p)!=1.
*/

/*
	Ricorda di inizializzare value ad ogni inserimento!
*/
struct base{
	base* next;
	unsigned long value;	
	//unsigned long num_of_element;
}

base* eratostene(unsigned long int n){
	/*
		Fa la cosa giusta
	*/
}


/*
	sarebbe meglio utilizzare una lista concatenata
	se esiste nella libreria del c
	stdlib.h
*/
/*
    Inizializzo all'esterno l'array con il numero giusto secondo l'approssimazione
    di elementi.

    La funzione ritorna il numero di elementi della base di fattori scremata
    dai non residui quadratici.

	Scorro tutti gli elementi della base di fattori 
	ritornata dalla funzione eratostene(), 
	quelli che hanno legendre(base.value)!=1 vanno eliminati
	Siamo sicuri che il primo elemento della base ritornata da eratostene sia
	2, quindi non dobbiamo eliminarlo.
*/
unsigned int factorBase(mpz_t N, unsigned int* fb_1, int dim1, unsigned int* fb_2, int dim2){
	
	int j=0;

	for(int i=0; i<dim1; i++){
		if(legendre(fb_1[i])==1){
			fb_2[j]=fb_1[i];
			++j;
		}
	}	

	return j;
}