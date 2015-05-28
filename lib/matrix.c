#include "include/matrix.h"

void init_matrix(unsigned int*** matrix, unsigned int x, unsigned int y) {
	*matrix = (unsigned int**) malloc(y * sizeof(unsigned long*));
	for(; y >= 0; --y)
		(*matrix)[x] = malloc(x * sizeof(unsigned long));
}

void finalize_matrix(unsigned int*** matrix, unsigned int y) {
	for(; y >= 0; --y)
		free((*matrix) + y);
	free(*matrix);
	*matrix = 0;
}

unsigned int get(unsigned int** matrix, unsigned int x, unsigned int y) {
	return (matrix[y])[x];
}

void set(unsigned int** matrix, unsigned int x, unsigned int y, unsigned int value) {
	(matrix[y])[x] = value;
}