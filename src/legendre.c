#include "../include/legendre.h"
#include <stdio.h>
#include <math.h>


// compile gcc legendre.c -std=c99 -lgmp -o legendre


/*
	All'interno del crivello quadratico, quando vado a scegliere la base di fattori, 
	devo eliminare quei primi p tali per cui il simbolo di legendre (N|p)!=1, 
	con N il numero da fattorizzare.
	Per ogni primo della base devo eseguire il controllo rispetto ad N.


	Sia p un numero primo ed a un numero intero, allora
	il simbolo di Legendre 
 
		    	{	1 se p|/a e a è un residuo quadratico modulo p
	(a|p) = |	0 se p|a
		    	{   -1 se p|/a e a non è un residuo quadratico modulo p

	valgono le proprietà seguenti:
	
	periodicità             (a|p) = (a mod p | p)  
	regola del 2            (2|p) = (-1)^((p^2-1)/8) 
	reciprocità quadratica  (p|q)(q|p) = (-1)^((p-1)*(q-1)/8)
	moltiplicatività        (ab|p) = (a|p)(b|p) 


	La base dovrebbe essere composta da interi relativamente "piccoli", nell'ordine di 
	exp((1/2+o(1)*sqrt(log(N)(loglog(N)))))
	per numeri come gli rsa non basta considerare un unsigned long long perché la quantità
	di numeri primi è nell'ordine di 10^5.

*/




void inizializza_pair(pair *p){
  p->sol1=0;
  p->sol2=0;
}



/*
  Costo O(p/2)
*/
int legendre(mpz_t n, unsigned int p ){


int sol;
mpz_t m;
mpz_init(m);
mpz_mod_ui(m,n,p); // m = n % p;

unsigned int mm = mpz_get_ui(m);

  for (unsigned j = 0; 2 * j < p ; ++j){
    if ((j*j) % p == mm) {
      return 1;
    }
  }

  return 0;
}




/*
  Calcola (n|p).
  s=sqrt(N)
  in solution metto le soluzioni dell'equazione x^2=n mod p
*/

/*
  Costo O(p)
*/
int legendre_sol(mpz_t n, unsigned int p, pair *sol ){
 
  // questa variabile contiene il risultato
  // 1 se n è residuo quadratico mod p, 0 altrimenti
  int leg = 0;
  /*
    Parte del prof per la ricerca delle soluzioni di 
     x^2=a mod p
  */
  mpz_t m; //     m = n % p;
  mpz_init(m);
  mpz_mod_ui(m,n,p);

  unsigned int mm = mpz_get_ui(m);
 
  // s = sup(sqrt(n))
  unsigned int n1 = mpz_get_ui(n);
  float s1 = sqrt(n1);
  s1 = ceil(s1);

  unsigned int s = (unsigned int)(s1);

 
  // dichiaro j fuori per poter riiniziare il secondo ciclo
  // dall'ultimo valore valutato dal primo
  unsigned j;

  /*
      Uso due cicli for per calcolare le due soluzioni dell'equazione x^2=n mod p
      Al primo ciclo mi fermo subito dopo aver trovato la prima soluzione.
      Nel secondo riparto dal valore sul quale mi ero precedentemente fermato.
      NB if (j!=0) serve solo per non entrare nel secondo ciclo se non si è passati dal primo
  */

    for (unsigned j = 0; 2 * j < p ; ++j){
      if ((j*j) % p == mm) {
        leg = 1;
      }
    }  

    for (j = 0;  j < p ; ++j){
      if ((j*j) % p == mm) {
        sol->sol1 = j;
        break;
      }
    }  

    if(j!=0){
      j++;
      for (j;  j < p ; ++j){
        if ((j*j) % p == mm) {
         sol->sol2 = j;
         break;
        }
      }
    }  

   return leg;
}





int main(){

    printf("ciao pippo1 \n");

	mpz_t n;
	mpz_init(n);
	mpz_set_ui(n,4);

    printf("ciao pippo2 \n");

    unsigned int p = 17;

    printf("ciao pippo3 \n");
    
    mpz_t s;
    mpz_init(s);
    mpz_sqrt(s,n);

    printf("ciao pippo4 \n");

    pair solution;

    printf("ciao pippo5 \n"); 

    inizializza_pair(&solution);

    printf("ciao pippo6 \n");

    unsigned int pippo=-3;

    printf("ciao pippo7 \n");

    pippo = legendre(n,p);

    printf("ciao pippo8 \n");

    printf("%d\n",pippo);

    pippo = legendre_sol(n,p,&solution);

    printf("ciao pippo9 \n");

    printf("%d\n",pippo);

    printf("%d\n",solution.sol1);
    printf("%d\n",solution.sol2);    

}
