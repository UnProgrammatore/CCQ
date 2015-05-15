/*
* QUESTO NON COMPILA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <gmp.h>

struct elem { // Very basic and non-reusable stack
	long long val;
	struct elem* next;
};

void add(struct elem** head, long long val) {
	struct elem* app = malloc(sizeof(struct elem));
	app->val = val;
	app->next = *head;
	*head = app;
}

long long pick(struct elem** head) {
	long long toret;
	struct elem* app;
	if(*head == NULL)
		return 0;
	else {
		toret = (*head)->val;
		app = *head;
		*head = (*head)->next;
		free(app);
		return toret;
	}
}

void master_procedure(int comm_size) {
	int i = 1;
	long long rec;
	int shit_happened;
	
	while(i < comm_size) {
		shit_happened = MPI_Recv(&rec, 1, MPI_LONG_LONG, i, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		
		if(shit_happened) {
			fprintf(stderr, "Recv failed");
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		if(rec == 0)
			++i;
		else
			printf("Factor: %lld\n", rec);
	}
}

void slave_procedure(int my_rank, int comm_size, mpz_t the_number) {
	long long from, to;
	long long to_send;
	int shit_happened;
	struct elem* head = NULL;

	mpz_t temp;
	mpzt_t from;
	mpz_t to;
	mpz_init(temp);
	mpz_init(from);
	mpz_init(to);

	mpz_root(temp, the_number, 2);
	mpz_div_ui(temp, temp, comm_size - 1);
	
	mpz_mul_ui(from, temp, my_rank - 1);
	mpz_mul_ui(from, temp, my_rank);	

	mpz_cmp_ui(from, 0) ? : mpz_set_ui(from, 1); // Because why not

	#pragma omp parallel shared(from, to)
	{
		int threads = omp_get_thread_num();
		int my_thread = omp_get_num_threads();
		
		mpz_t from_thread;
		mpz_t to_thread;
		mpz_t module;
		mpz_init(from_thread);
		mpz_init(to_thread);
		mpz_init(module);

		mpz_sub(to_thread, to, from);
		mpz_set(from_thread, to_thread);

		mpz_div_ui(to_thread, to_thread, threads); // CONTROLLA E FINISCI DA QUI
		mpz_mul_ui(to_thread, to_thread, my_thread + 1);

		mpz_div_ui(from_thread, from_thread, threads);
		mpz_mul_ui(from_thread, from_thread, my_thread);

		mpz_add(from_thread, from_thread, from);
		mpz_add(to_thread, to_thread, from);

		while(mpz_cmp(from_thread, to_thread) <= 0) {
			// DOPO QUI NON HO ANCORA FINITO
			if(the_number % (from + i) == 0) {
				#pragma omp critical
				{
					add(&head, from + i);
					add(&head, the_number / (from + i));
				}
				//QUISOMMA
			}
		}
	}
	do {
		to_send = pick(&head);
		shit_happened = MPI_Ssend(&to_send, 1, MPI_LONG_LONG, 0, 0, MPI_COMM_WORLD);

		if(shit_happened) {
			fprintf(stderr, "Send failed");
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

	}while(to_send != 0);
}

int main(int argc, char** argv) {
	int my_rank, comm_size;
	mpz_t the_number;

	mpz_init(the_number);

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	if(argc <= 1) {
		fprintf(stderr, "Missing number as argument");
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	else
		mpz_set_str(the_number, argv[1], 10); // 10 is the base

	if(my_rank == 0)
		master_procedure(comm_size);
	else
		slave_procedure(my_rank, comm_size, the_number);

	MPI_Finalize();
	return 0;
}