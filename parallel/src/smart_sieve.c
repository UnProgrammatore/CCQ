#include "../include/smart_sieve.h"
#include "../include/linear_algebra.h"
#include "../include/matrix.h"

unsigned int smart_sieve(mpz_t n,
			 unsigned int* factor_base,
			 unsigned int base_dim,
			 pair* solutions,
			 mpz_t begin,
			 unsigned int interval,
			 unsigned int block_size,
			 unsigned int max_fact) {
  mpz_t end;
  mpz_init(end);

  // questo processo mpz deve prendere A in [begin, end]
  mpz_add_ui(end, begin, interval); // end = begin + interval

  mpz_t n_root;
  mpz_init(n_root);
  mpz_sqrt(n_root, n);

  int stop_flag = 0;
  char stop_signal;
  MPI_Request request;
  MPI_Status status;
  MPI_Irecv(&stop_signal, 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &request);
  
  /* Inizio della parte di codice eseguita da ogni thread */
  #pragma omp parallel
  { 

    /* Dichiarazione dei buffer per la trasmissione */
    unsigned int* buffer; // buffer per le fattorizzazioni
    init_vector(&buffer, base_dim); // dim+1 per il flag di controllo della fine
    unsigned char* buffer_as; // buffer per (A + s)
    buffer_as = malloc(sizeof(unsigned char) * BUFFER_DIM);

    /* Dichiarazione strutture dati per raccolta risultati */
    unsigned int** exponents; // Matrice temporanea degli esponenti
    init_matrix(&exponents, block_size, base_dim);
    word** used_rows; // Vettore che segna quali esponenti sono stati inizializzati
    init_matrix_l(&used_rows, 1, (block_size / N_BITS) + 1);
    mpz_t* evaluated_poly; // Vettore temporaneo dei Q(A) valutati
    init_vector_mpz(&evaluated_poly, block_size);
    mpz_t* As; // A + s temporanei
    init_vector_mpz(&As, block_size);

    // Ogni thread calcola per A in [begin_thread, end_thread]
    mpz_t begin_thread;
    mpz_init(begin_thread);
    mpz_t end_thread;
    mpz_init(end_thread);

    mpz_t last_block;
    mpz_init(last_block);
    mpz_t end_block;
    mpz_init(end_block);

    mpz_t intermed; // Valore appoggio
    mpz_init(intermed);

    mpz_t A; // Valore di A di Q(A) = (A + s)^2
    mpz_init(A);

    mpz_t l;
    mpz_init(l);

    mpz_t j;
    mpz_init(j);

    mpz_t begin_solution1;
    mpz_init(begin_solution1);
    mpz_t begin_solution2;
    mpz_init(begin_solution2);

    // Indice per accedere a Q(A) memorizzato in posizione (A - offset) tale che index < block_size
    mpz_t index_mpz; 
    mpz_init(index_mpz);
    unsigned int index;
   
    int n_bytes; // num bit per vettore righe usate nella matrice

    unsigned int i; // Indice generico
    unsigned char go_on = 1;
    unsigned int h; // Usato per copiare gli base_dim esponenti

    unsigned int fact_count = 0; // Numero fattorizzazioni trovate
    unsigned long k; // Indici generici

    /*******************************************************************************************************/

    max_fact += base_dim; // Le k + n fattorizzazioni da trovare

    int threads = omp_get_num_threads();
    int thread_id = omp_get_thread_num();

    unsigned int dom_decomp = interval / (threads);
    mpz_add_ui(begin_thread, begin, dom_decomp * thread_id); // begin_thread = begin + (dom_decomp * thread_id)
    mpz_add_ui(end_thread, begin_thread, dom_decomp); // end_thread = begin_thread + dom_decomp
    //gmp_printf("begin=%Zd, end=%Zd, interval=%d, block_size=%d\n", begin, end, interval, block_size);
    //gmp_printf("%d) begin_thread=%Zd, end_thread=%Zd, dom_decmp=%d\n", thread_id, begin_thread, end_thread, dom_decomp);
    //printf("###originali\n");
    mpz_pair * solutions_ = malloc(sizeof(mpz_pair) * base_dim);
    for(i = 0; i < base_dim; ++i) {
      mpz_init(solutions_[i].sol1);
      mpz_init(solutions_[i].sol2);

      mpz_set_ui(solutions_[i].sol1, solutions[i].sol1);
      mpz_set_ui(solutions_[i].sol2, solutions[i].sol2);

      //gmp_printf("x_%d = %Zd, ", factor_base[i], solutions_[i].sol1);
      //gmp_printf("y_%d = %Zd\n", factor_base[i], solutions_[i].sol2);
      
    }
    //printf("###ricalcolate\n"); 
    if(mpz_cmp_ui(begin_thread, 0) != 0)
      for(i = 0; i < base_dim; ++i) {
	//unsigned int old1 = solutions_[i].sol1;
	//unsigned int old2 = solutions_[i].sol2;
	
	//unsigned int f = (startfrom - solutions[i].sol1) / factor_base[i] + 1;
	//solutions[i].sol1 = solutions[i].sol1 + f * factor_base[i];
	//f = (startfrom - solutions[i].sol2) / factor_base[i] + 1;
	//solutions[i].sol2 = solutions[i].sol2 + f * factor_base[i];
	
	while(mpz_cmp(solutions_[i].sol1, begin_thread) < 0)
	  mpz_add_ui(solutions_[i].sol1, solutions_[i].sol1, factor_base[i]);
	while(mpz_cmp(solutions_[i].sol2, begin_thread) < 0)
	  mpz_add_ui(solutions_[i].sol2, solutions_[i].sol2, factor_base[i]);

	//gmp_printf("x_%d = %Zd, ", factor_base[i], solutions_[i].sol1);
	//gmp_printf("y_%d = %Zd\n", factor_base[i], solutions_[i].sol2);
      }
    //printf("###fine calcolo soluzioni \n");

    mpz_sub_ui(last_block, end_thread, block_size); // last_block = end_thread - block_size

    // for(l = begin_thread; l < last_block && go_on; l += block_size)
    for(mpz_set(l, begin_thread); (mpz_cmp(l, last_block) < 0) && go_on && (stop_flag == 0); mpz_add_ui(l, l, block_size)) {
      for(i = 0; i < ((block_size / N_BITS) + 1); ++i) { // Reset righe usate
	set_matrix_l(used_rows, 0, i, 0);
      }

      mpz_add_ui(end_block, l, block_size); // end_block = l + block_size
      //gmp_printf("l=%Zd < %Zd [%Zd, %Zd)\n", l, last_block, l, end_block);

      for(i = 0; i < block_size; ++i) { // Calcolo Q(A) e (A + s) per A in [l, l + block_size]
	mpz_add_ui(A, l, i); // A = i + l
	mpz_add(intermed, n_root, A); // A + s
	mpz_set(As[i], intermed);
	mpz_mul(intermed, intermed, intermed); // (A + s)^2
	mpz_sub(evaluated_poly[i], intermed, n);

	//gmp_printf("Q(%Zd)=%Zd, ", A, evaluated_poly[i]);
      }
      //printf("\n");

      for(i = 0; i < base_dim && go_on; ++i) {

	/* Sieve con Xp */
	// for(j = solutions_[i].sol1; j < end_block && go_on; j += factor_base[i])
	for(mpz_set(j, solutions_[i].sol1); (mpz_cmp(j, end_block) < 0) && go_on; mpz_add_ui(j, j, factor_base[i])) {
          //gmp_printf("\txp) j=%Zd < %Zd [+=%d](j, end_block)\n", j, end_block, factor_base[i]);
	  mpz_sub(index_mpz, j, l);
	  index = mpz_get_ui(index_mpz); // Siccome (j - l) < block_size è un uint sicuro

	  while(mpz_divisible_ui_p(evaluated_poly[index], factor_base[i])) {
            //gmp_printf("\t\tQ(A) = %Zd / %d = ", evaluated_poly[index], factor_base[i]);
	    if(get_k_i(used_rows, 0, index) == 0) { // Se non sono mai stati usati gli esponenti
	      for(k = 0; k < base_dim; ++k)
		set_matrix(exponents, index, k, 0);
	      set_k_i(used_rows, 0, index, 1);
	    }

	    set_matrix(exponents, index, i, get_matrix(exponents, index, i) + 1); // ++exponents[j][i]
	    mpz_divexact_ui(evaluated_poly[index], evaluated_poly[index], factor_base[i]); // Q(A) = Q(A) / p
	    //gmp_printf("%Zd  (poly[%d])\n", evaluated_poly[index], index);
	  }
	}
	mpz_set(solutions_[i].sol1, j); // solutions[i].sol1 = j; // Al prossimo giro ricominciamo da dove abbiamo finito

	/* Sieve con Yp */
	// for(j = solutions_[i].sol2; j < end_block && go_on; j += factor_base[i])
	for(mpz_set(j, solutions_[i].sol2); factor_base[i] != 2 && (mpz_cmp(j, end_block) < 0) && go_on; mpz_add_ui(j, j, factor_base[i])) {
          //gmp_printf("\txp) j=%Zd < %Zd [+=%d](j, end_block)\n", j, end_block, factor_base[i]);
	  mpz_sub(index_mpz, j, l);
	  index = mpz_get_ui(index_mpz); // Siccome (j - l) < block_size è un uint sicuro

	  while(mpz_divisible_ui_p(evaluated_poly[index], factor_base[i])) {
            //gmp_printf("\t\tQ(A) = %Zd / %d = ", evaluated_poly[index], factor_base[i]);
	    if(get_k_i(used_rows, 0, index) == 0) { // Se non sono mai stati usati gli esponenti
	      for(k = 0; k < base_dim; ++k)
		set_matrix(exponents, index, k, 0);
	      set_k_i(used_rows, 0, index, 1);
	    }

	    set_matrix(exponents, index, i, get_matrix(exponents, index, i) + 1); // ++exponents[j][i]
	    mpz_divexact_ui(evaluated_poly[index], evaluated_poly[index], factor_base[i]); // Q(A) = Q(A) / p
	    //gmp_printf("%Zd  (poly[%d])\n", evaluated_poly[index], index);
	  }
	}
	mpz_set(solutions_[i].sol2, j); // solutions[i].sol2 = j; // Al prossimo giro ricominciamo da dove abbiamo finito
      }

      // Spedisco le fattorizzazioni trovate in questo blocco
      for(i = 0; i < block_size; ++i) {
	if(mpz_cmp_ui(evaluated_poly[i], 1) == 0) {
	  ++fact_count;
	  for(k = 0; k < base_dim; ++k)
	    buffer[k] = get_matrix(exponents, i, k);

	  /* MPI_Send */
          #pragma omp critical 
	  {
	    MPI_Send(buffer, base_dim, MPI_UNSIGNED, 0, ROW_TAG, MPI_COMM_WORLD);
	    n_bytes = (mpz_sizeinbase(As[i], 2) + 7) / 8;
	    *buffer_as = 0;
	    mpz_export(buffer_as, NULL, 1, 1, 1, 0, As[i]);
	    MPI_Send(buffer_as, n_bytes, MPI_UNSIGNED_CHAR, 0, AS_TAG, MPI_COMM_WORLD);

	    if(stop_flag == 0)
	      MPI_Test(&request, &stop_flag, &status);
	  }
	}
      }
    }	      
  }
  return stop_flag;
}
