
#include<gmp.h>
#include<stdlib.h>
#include<math.h>
#include<stdio.h>

// compile gcc Legendre.c -std=c99 -lgmp -o Legendre


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




typedef
struct pair_{
	unsigned int sol1;
	unsigned int sol2;
} pair;

void inizializza_pair(pair p){
	p.sol1=0;
	p.sol2=0;
}

/*
	Calcola (n|p).
	s=sqrt(N)
	in solution metto le soluzioni dell'equazione x^2=n mod p
*/
unsigned int legendre(mpz_t n, unsigned int p, pair sol ){


    int res=0;

    mpz_t m; 
    mpz_init(m);
    mpz_mod_ui(m,n,p);  // m = n%p ;
    // Data la simmetria delle soluzioni di "x^2 = n mod p" quando "p"
    // è dispari, è sufficiente esaminare i valori di "x"
    // nell'intervallo [0, (p - 1) / 2]
    for (unsigned j = 0 ; 2*j<p ; ++j)
      if (mpz_cmp_ui(m,(j*j) % p)==0) {  // se j*j%p == m
	// Se si entra qui vuol dire che "n" è un residuo quadratico
	// modulo "p", e quindi "p" fa parte della "base di fattori".
	// Le due soluzioni di "Q(A) = n mod p" si possono calcolare a
	// partire dalle soluzioni di "x^2 = n mod p".  Bisogna
	// accertarsi che la soluzione calcolata cada nell'intervallo
	// [0, p - 1], e per questo sono necessarie le operazioni di
	// modulo.
    
      res = 1;

      /*
			converto  p, j in mpz per trovare sol1 e sol2.
			 	  sol.sol1 = (p - ((s - j) % p)) % p; 
				  sol.sol2 = (p - ((j + s) % p)) % p;

      */

      mpz_t pm;
      mpz_init(pm);

      mpz_t jm;
      mpz_init(jm);

      mpz_set_ui(pm,p);
      mpz_set_ui(jm,j);

      mpz_t sol1m;
      mpz_init(sol1m);

	  mpz_t sol2m;
      mpz_init(sol2m);

      mpz_t smj, spj,p1,p2; // s-j   s+j
      mpz_init(smj);
      mpz_init(spj);
      mpz_init(p1);
      mpz_init(p2);

      mpz_sub(smj,s,j);
      mpz_add(spj,s,j);

      mpz_mod(smj,smj,pm);   // smj = smj mod p
      mpz_mod(spj,spj,pm);	 //	spj = spj mod p

      mpz_sub(p1,pm,smj);    // p1 = p-(s-j)
      mpz_sub(p2,pm,spj);    // p2 = p-(s+j)

      mpz_mod(p1,p1,pm);	 // p1 = p1 mod p
      mpz_mod(p1,p1,pm);     // p2 = p2 mod p

      sol.sol1=mpz_get_ui(p1);
      sol.sol2=mpz_get_ui(p2);
	}
    
    return res;

}


int main(){

    printf("ciao pippo1 \n");

	mpz_t n;
	mpz_init(n);
	mpz_set_ui(n,1234);

    printf("ciao pippo2 \n");

	mpz_t p;
	mpz_init(p);
	mpz_set_ui(p,23);

    printf("ciao pippo3 \n");
    
    mpz_t s;
    mpz_init(s);
    mpz_sqrt(s,n);

    printf("ciao pippo4 \n");

    pair solution;

    printf("ciao pippo5 \n"); 

    inizializza_pair(solution);

    printf("ciao pippo6 \n");

    unsigned int pippo;

    printf("ciao pippo7 \n");

	pippo = legendre(n,p,s,solution);

    printf("ciao pippo8 \n");

    printf("%d\n",pippo);


}
