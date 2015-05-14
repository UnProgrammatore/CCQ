#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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

void slave_procedure(int my_rank, int comm_size, long long the_number) {
	long long from, to;
	long long to_send;
	int shit_happened;
	struct elem* head = NULL;

	from = ((((long long int) sqrt((double) the_number)) / (comm_size)) + 1) * (my_rank - 1);
	to = ((((long long int) sqrt((double) the_number)) / (comm_size)) + 1) * (my_rank); // TODO: better square root

	from = from == 0 ? 1 : from; // Because why not

	for(; from < to; ++from) {
		if(the_number % from == 0) {
			add(&head, from);
			add(&head, the_number / from);
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
	long long the_number;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	if(argc <= 1) {
		fprintf(stderr, "Missing number as argument");
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	else
		the_number = atoll(argv[1]);

	if(my_rank == 0)
		master_procedure(comm_size);
	else
		slave_procedure(my_rank, comm_size, the_number);

	MPI_Finalize();
	return 0;
}