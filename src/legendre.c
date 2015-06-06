#include "legendre.h"


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




void inizializza_pair(pair p){
  p.sol1=0;
  p.sol2=0;
}



int legendre(mpz_t n, unsigned int p ){


int sol;
mpz_t m;
mpz_init(m)
mpz_mod_ui(m,n,p); // m = n % p;

  for (unsigned j = 0; 2 * j < p ; ++j){
    if ((j*j) % p == m) {
      return 1;
    }
  }

}

/*
int legendre(mpz_t n, unsigned int p ){

  p = (p-1)/2;   // sono sicuro sia intero
  mpz_t res;
  mpz_init(res);
  mpz_pow_ui(res, n, p);  // res = a^((p-1)/2)

  mpz_mod_ui(res, res, p);

  if(mpz_cmp_ui(res,0))
    return 0;

  if(mpz_cmp_ui(res,1))
    return 1;

  return -1;

}
*/



/*
  Calcola (n|p).
  s=sqrt(N)
  in solution metto le soluzioni dell'equazione x^2=n mod p
*/


int legendre_sol(mpz_t n, unsigned int p, pair sol ){
 
  int leg = 0;
  /*
    Parte del prof per la ricerca delle soluzioni di 
     x^2=a mod p
  */
  mpz_t m; //     m = n % p;
  mpz_init(m);
  mpz_mod_ui(m,n,p);

  mpz_t s;
  mpz_init(s);
  mpz_sqrt(s,n);


    // Data la simmetria delle soluzioni di "x^2 = n mod p" quando "p"
    // è dispari, è sufficiente esaminare i valori di "x"
    // nell'intervallo [0, (p - 1) / 2]
    for (unsigned j = 0; 2 * j < p ; ++j){
      if ((j*j) % p == m) {
      leg = 1;  

      // Se si entra qui vuol dire che "n" è un residuo quadratico
      // modulo "p", e quindi "p" fa parte della "base di fattori".
      // Le due soluzioni di "Q(A) = n mod p" si possono calcolare a
      // partire dalle soluzioni di "x^2 = n mod p".  Bisogna
      // accertarsi che la soluzione calcolata cada nell'intervallo
      // [0, p - 1], e per questo sono necessarie le operazioni di
      // modulo.
      sol.sol1 = (p - ((s - j) % p)) % p;
      sol.sol2 = (p - ((j + s) % p)) % p;

      return leg;
   }
}


}

/*3
int main(){

    printf("ciao pippo1 \n");

	mpz_t n;
	mpz_init(n);
	mpz_set_ui(n,22);

    printf("ciao pippo2 \n");

    unsigned int p = 23;

    printf("ciao pippo3 \n");
    
    mpz_t s;
    mpz_init(s);
    mpz_sqrt(s,n);

    printf("ciao pippo4 \n");

    pair solution;

    printf("ciao pippo5 \n"); 

    //inizializza_pair(solution);

    printf("ciao pippo6 \n");

    unsigned int pippo;

    printf("ciao pippo7 \n");

    pippo = legendre(n,p);

    printf("ciao pippo8 \n");

    printf("%d\n",pippo);


}
*/