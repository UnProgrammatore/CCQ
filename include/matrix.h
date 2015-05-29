#ifndef MATRIX_GUARD
#define MATRIX_GUARD

void init_matrix(unsigned int*** matrix, unsigned int x, unsigned int y);
void finalize_matrix(unsigned int*** matrix, unsigned int x);
unsigned int get_matrix(unsigned int** matrix, unsigned int x, unsigned int y);
void set_matrix(unsigned int** matrix, unsigned int x, unsigned int y, unsigned int value);

void init_matrix_long(unsigned long*** matrix, unsigned int x, unsigned int y);
void finalize_matrix_long(unsigned long*** matrix, unsigned int x);
unsigned int get_matrix_long(unsigned long** matrix, unsigned int x, unsigned int y);
void set_matrix_long(unsigned long** matrix, unsigned int x, unsigned int y, unsigned int value);

#endif // MATRIX_QUARD