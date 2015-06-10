#include <gmp.h>
#include "../include/matrix.h"
#include <stdio.h>

int main() {
	mpz_t** matrice;
	init_matrix_mpz(&matrice, 5, 5);
	unsigned int ui;
	mpz_t mp;
	mpz_init(mp);
	mpz_set_str(mp, "258584478548548754825896954865781478158", 10);
	printf("Sto per settare l'mpz\n");
	set_matrix_mpz(matrice, 2, 2, mp);

	mpz_t val;
	mpz_init(val);
	printf("Sto per leggere l'mpz\n");
	get_matrix_mpz(val, matrice, 2, 2);
	gmp_printf("In 2, 1 ho %Zd\n", val);
	
	printf("Sto per scrivere un altro mpz\n");
	mpz_set_str(mp, "456", 10);
	set_matrix_mpz(matrice, 3, 2, mp);

	printf("Sto per leggere l'altro mpz\n");
	get_matrix_mpz(val, matrice, 3, 2);
	gmp_printf("Ho letto %Zd\n", val);

	printf("Scrivo sul bordo\n");
	mpz_set_str(mp, "111", 10);
	set_matrix_mpz(matrice, 4, 4, mp);

	printf("Leggo dal bordo\n");
	get_matrix_mpz(val, matrice, 4, 4);
	gmp_printf("Letto dal bordo %Zd\n", val);

	/*
	ui = 5;

	printf("Coso ui altro\n");
	mpz_set_ui((matrice[2])[1], ui);

	printf("Sto per settare l'ui\n");
	set_matrix_mpz_ui(matrice, 2, 1, ui);
	
	printf("Sto per leggere l'ui\n");
	get_matrix_mpz(val, matrice, 2, 1);
	gmp_printf("In 2, 2 ho %Zd\n", val);*/

}