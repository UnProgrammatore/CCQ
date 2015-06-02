#include "matrix.h"
#include <stdlib.h>

#include <stdio.h>

void init_matrix(unsigned int*** matrix, unsigned int x, unsigned int y) {
	*matrix = (unsigned int**) malloc(y * sizeof(unsigned int*));
	for(; y > 0; --y) {
		(*matrix)[y - 1] = malloc(x * sizeof(unsigned int));
	}
}

void finalize_matrix(unsigned int*** matrix, unsigned int y) {
	for(; y > 0; --y) {
		free((*matrix)[y - 1]);
	}
	free(*matrix);
	*matrix = 0;
}

unsigned int get_matrix(unsigned int** matrix, unsigned int x, unsigned int y) {
	return (matrix[y])[x];
}

void set_matrix(unsigned int** matrix, unsigned int x, unsigned int y, unsigned int value) {
	(matrix[y])[x] = value;
}


void init_matrix_l(unsigned long*** matrix, unsigned int x, unsigned int y) {
	*matrix = (unsigned long**) malloc(y * sizeof(unsigned long*));
	for(; y > 0; --y) {
		(*matrix)[y - 1] = malloc(x * sizeof(unsigned long));
	}
}

void finalize_matrix_l(unsigned long*** matrix, unsigned int y) {
	for(; y > 0; --y) {
		free((*matrix)[y - 1]);
	}
	free(*matrix);
	*matrix = 0;
}

unsigned long get_matrix_l(unsigned long** matrix, unsigned int x, unsigned int y) {
	return (matrix[y])[x];
}

void set_matrix_l(unsigned long** matrix, unsigned int x, unsigned int y, unsigned int value) {
	(matrix[y])[x] = value;
}
