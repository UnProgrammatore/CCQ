#include <stdio.h>
#include <gmp.h>

typedef struct sol_pair_ {
	long unsigned sol1;
	long unsigned sol2;
} sol_pair;

short uguali_modulo_mpz_ui(long unsigned f1, mpz_t f2, long unsigned modulo) {
	mpz_t appoggio;
	mpz_init(appoggio);

	mpz_mod_ui(appoggio, f2, modulo);

	return (f1 % modulo) == mpz_get_ui(appoggio);
}

sol_pair calcola_cosi(mpz_t radice, long unsigned j, long unsigned p) {
	mpz_t appoggio;
	mpz_init(appoggio);

	sol_pair pair;

	mpz_sub_ui(appoggio, radice, j); // appoggio = radice - j;
	mpz_mod_ui(appoggio, appoggio, p);

	pair.sol1 = mpz_get_ui(appoggio);

	mpz_add_ui(appoggio, radice, j);
	mpz_mod_ui(appoggio, appoggio, p);

	pair.sol2 = mpz_get_ui(appoggio);

	return pair;

}

long unsigned base_fattori(mpz_t numero, mpz_t radice, long unsigned base_fattori[], sol_pair soluzioni[], long unsigned primi[], long unsigned quanti_primi) {
	base_fattori[0] = 2;
	mpz_t soluzione;
	mpz_init(soluzione);
	
	mpz_add_ui(soluzione, radice, 1);
	mpz_mod_ui(soluzione, soluzione, 2);

	soluzioni[0].sol1 = soluzioni[0].sol2 = mpz_get_ui(soluzione);


	long unsigned k = 1;
	long unsigned i = 1;
	long unsigned j;

	sol_pair pair;

	for(i = 1; i < quanti_primi; ++i) {
		
		for(j = 0; 2 * j < primi[i]; ++j) {
			
			if(uguali_modulo_mpz_ui(j*j, numero, primi[i])) {
				
				pair = calcola_cosi(radice, j, primi[i]);

				soluzioni[k].sol1 = (primi[i] - pair.sol1) % primi[i];
				soluzioni[k].sol2 = (primi[i] - pair.sol2) % primi[i];

				base_fattori[k] = primi[i];

				++k;
			}

		}
	}
}