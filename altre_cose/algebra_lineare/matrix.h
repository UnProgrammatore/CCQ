#ifndef MATRIX_GUARD
#define MATRIX_GUARD

void init_matrix(unsigned int*** matrix, unsigned int x, unsigned int y);
void finalize_matrix(unsigned int*** matrix, unsigned int x);
unsigned int get_matrix(unsigned int** matrix, unsigned int x, unsigned int y);
void set_matrix(unsigned int** matrix, unsigned int x, unsigned int y, unsigned int value);

void init_matrix_l(unsigned long*** matrix, unsigned int x, unsigned int y);
void finalize_matrix_l(unsigned long*** matrix, unsigned int x);
unsigned long get_matrix_l(unsigned long** matrix, unsigned int x, unsigned int y);
void set_matrix_l(unsigned long** matrix, unsigned int x, unsigned int y, unsigned long value);

#endif // MATRIX_QUARD
