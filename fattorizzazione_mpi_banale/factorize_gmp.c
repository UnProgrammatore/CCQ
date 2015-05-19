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
	mpz_t val;
	struct elem* next;
};

void add(struct elem** head, mpz_t val) {
	struct elem* app = malloc(sizeof(struct elem));
	mpz_init(app->val);
	mpz_set(app->val, val); // app->val = val;
	app->next = *head;
	*head = app;
}

mpz_t pick(struct elem** head) {
	mpz_t toret;
	mpz_init(toret);

	struct elem* app;
	
	if(*head == NULL)
		mpz_set_ui(toret, 0); // toret = 0;
	else {
		mpz_set(toret, (*head)->val); // toret = (*head)->val;
		app = *head;
		*head = (*head)->next;
		mpz_finalize(app->val);
		free(app);
	}
	return toret;
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
	int shit_happened;
	struct elem* head = NULL;

	mpz_t temp;
	mpzt_t from;
	mpz_t to;
	mpz_init(temp);
	mpz_init(from);
	mpz_init(to);

	mpz_root(temp, the_number, 2); // temp = sqrt(the_number);
	mpz_div_ui(temp, temp, comm_size - 1); // temp = temp / (comm_size - 1);
	
	mpz_mul_ui(from, temp, my_rank - 1); // from = temp * (my_rank - 1);
	mpz_mul_ui(to, temp, my_rank); // to = temp * my_rank;

	mpz_cmp_ui(from, 0) ? : mpz_set_ui(from, 1); // from == 0 ? from = 1 : ;

	#pragma omp parallel shared(from, to)
	{
		int threads = omp_get_thread_num();
		int my_thread = omp_get_num_threads();
		
		mpz_t from_thread;
		mpz_t to_thread;
		mpz_t divided;
		mpz_t to_send;
		mpz_init(from_thread);
		mpz_init(to_thread);
		mpz_init(divided);
		mpz_init(to_send);

		mpz_sub(to_thread, to, from); // to_thread = to - from;
		mpz_set(from_thread, to_thread); // from_thread = to_thread;

		mpz_div_ui(to_thread, to_thread, threads); // to_thread = to_thread / threads;
		mpz_mul_ui(to_thread, to_thread, my_thread + 1); // to_thread = to_thread * (my_thread + 1);

		mpz_div_ui(from_thread, from_thread, threads); // from_thread = from_thread / threads;
		mpz_mul_ui(from_thread, from_thread, my_thread); // from_thread = from_thread * my_thread;

		mpz_add(from_thread, from_thread, from); // from_thread = from_thread + from;
		mpz_add(to_thread, to_thread, from); // to_thread = to_thread + from;

		while(mpz_cmp(from_thread, to_thread) <= 0) {
			
			if(mpz_divisible_p(the_number, from_thread)) {

				mpz_divexact(divided, the_number, from_thread); // divided = the_number / from_thread; // Only works if the_number % from_thread == 0;
				
				#pragma omp critical
				{
					add(&head, from_thread);
					add(&head, divided);
				}
				mpz_add_ui(from_thread, from_thread, 1); // ++from_thread;
			}
		}
	}

	// TODO IMPORTANT: make work with gmp
	do {
		to_send = pick(&head);
		int how_many_bytes = mpz_sizeinbase(to_send, 2) + 
		
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