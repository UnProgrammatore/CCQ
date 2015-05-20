#ifndef MATH_GMP_GUARD
#define MATH_GMP_GUARD

#include <gmp.h>

void mcd(mpz_t output, mpz_t in1, mpz_t in2);
void mul_mod_n(mpz_t output, mpz_t mul1, mpz_t mul2, mpz_t mod);
void pow_mod_n(mpz_t output, mpz_t base, mpz_t exp, mpz_t mod);

#endif MATH_GMP_GUARD