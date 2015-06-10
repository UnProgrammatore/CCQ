#ifndef MATRIX_GUARD
#define MATRIX_GUARD

#include <gmp.h>

void init_matrix(unsigned int*** matrix, unsigned int x, unsigned int y);
void finalize_matrix(unsigned int*** matrix, unsigned int x);
unsigned int get_matrix(unsigned int** matrix, unsigned int x, unsigned int y);
void set_matrix(unsigned int** matrix, unsigned int x, unsigned int y, unsigned int value);

void init_matrix_l(unsigned long*** matrix, unsigned int x, unsigned int y);
void finalize_matrix_l(unsigned long*** matrix, unsigned int x);
unsigned long get_matrix_l(unsigned long** matrix, unsigned int x, unsigned int y);
void set_matrix_l(unsigned long** matrix, unsigned int x, unsigned int y, unsigned long value);

void init_matrix_mpz(mpz_t*** matrix, unsigned int x, unsigned int y);
void get_matrix_mpz(mpz_t value, mpz_t** matrix, unsigned int x, unsigned int y);
void set_matrix_mpz(mpz_t** matrix, unsigned int x, unsigned int y, mpz_t value);

#endif // MATRIX_QUARD
