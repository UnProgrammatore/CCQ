#include "../include/base_fattori.h"

short uguali_modulo_mpz_ui(long unsigned f1, mpz_t f2, long unsigned modulo) {
	mpz_t appoggio;
	mpz_init(appoggio);
	unsigned long appoggio_ui;

	mpz_mod_ui(appoggio, f2, modulo);
	appoggio_ui = mpz_get_ui(appoggio);
	
	mpz_clear(appoggio);	
	
	return (f1 % modulo) == appoggio_ui;
}

pair calcola_soluzioni(mpz_t radice, long unsigned j, long unsigned p) {
	mpz_t appoggio;
	mpz_init(appoggio);

	pair pair;

	mpz_sub_ui(appoggio, radice, j); // appoggio = radice - j;
	mpz_mod_ui(appoggio, appoggio, p);

	pair.sol1 = mpz_get_ui(appoggio);
	pair.sol1 = (p - pair.sol1) % p;

	mpz_add_ui(appoggio, radice, j);
	mpz_mod_ui(appoggio, appoggio, p);

	pair.sol2 = mpz_get_ui(appoggio);
	pair.sol2 = (p - pair.sol2) % p;

	return pair;

}

long unsigned base_fattori(mpz_t numero, mpz_t radice, unsigned int base_fattori[], pair soluzioni[], 
			   unsigned int primi[], long unsigned quanti_primi) {
	base_fattori[0] = 2;
	mpz_t soluzione;
	mpz_init(soluzione);
	
	mpz_add_ui(soluzione, radice, 1);
	mpz_mod_ui(soluzione, soluzione, 2);

	soluzioni[0].sol1 = soluzioni[0].sol2 = mpz_get_ui(soluzione);


	long unsigned k = 1;
	long unsigned i = 1;
	long unsigned j;

	pair pair;

	for(i = 1; i < quanti_primi; ++i) {
		
		for(j = 0; 2 * j < primi[i]; ++j) {
			
			if(uguali_modulo_mpz_ui(j*j, numero, primi[i])) {
				
				soluzioni[k] = calcola_soluzioni(radice, j, primi[i]);
				base_fattori[k] = primi[i];

				++k;
			}

		}
	}

	//mpz_clear(soluzione);

	return k;
}
