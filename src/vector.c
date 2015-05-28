#include "../include/vector.h"

void init_vector(unsigned int** vector, unsigned int n) {
	*vector = (unsigned int*) malloc(sizeof(unsigned int) * n);
}


void finalize_vector(unsigned int** vector) {
	free(*vector);
	*vector = 0;
}