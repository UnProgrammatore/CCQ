#ifndef MATRIX_GUARD
#define MATRIX_GUARD

void init_matrix(unsigned int*** matrix, unsigned int x, unsigned int y);
void finalize_matrix(unsigned int*** matrix, unsigned int x);
unsigned int get(unsigned int** matrix, unsigned int x, unsigned int y);
void set(unsigned int** matrix, unsigned int x, unsigned int y, unsigned int value);

#endif // MATRIX_QUARD