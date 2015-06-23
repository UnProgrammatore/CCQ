#include "../include/smart_sieve.h"
#include "../include/linear_algebra.h"
#include "../include/matrix.h"

#include <omp.h>

unsigned int smart_sieve(
			 mpz_t n,
			 unsigned int* factor_base,
			 unsigned int base_dim,
			 pair* sol,
			 unsigned int M,
			 unsigned int max_fact,
			 unsigned int intervals,
			 unsigned int start
			 ) {

    #pragma omp parallel
    {
      
      //int flag = 0;
      
      int threads = omp_get_num_threads();
      int thread_id = omp_get_thread_num();

      //printf("thread %d/%d) s=%d M=%d\n",  thread_id, threads, start, M);

      unsigned int dom_decomp = (M-start) / (threads);
      unsigned int startfrom = start + (dom_decomp * thread_id);
      unsigned int poly_val_num = startfrom + dom_decomp;
     

      //printf("%d/%d) s=%d M=%d [%d, %d] i=%d intervals=%d dom_dec=%d\n", thread_id, threads, start, M, startfrom, poly_val_num, dom_decomp, intervals, dom_decomp);

      int how_many_bytes;

      unsigned int* buffer;
      init_vector(&buffer, base_dim+1);
      unsigned char* buffer_as;
      buffer_as = malloc(sizeof(unsigned char) * BUFFER_DIM);

      unsigned int i; // Indice generico

      unsigned char go_on = 1;

      unsigned int h; // Usato per copiare gli base_dim esponenti

      unsigned int** expo2; // Matrice temporanea degli esponenti
      init_matrix(&expo2, intervals, base_dim);
      word** is_used_expo2; // Vettore che segna quali esponenti sono stati inizializzati
      init_matrix_l(&is_used_expo2, 1, (intervals / N_BITS) + 1);

      mpz_t n_root;
      mpz_t intermed; // Valore appoggio
      mpz_init(n_root);
      mpz_init(intermed);
      mpz_sqrt(n_root, n);

      unsigned int fact_count = 0; // Numero fattorizzazioni trovate
      unsigned int j, k, l; // Indici generici
      mpz_t* evaluated_poly; // Vettore temporaneo dei Q(A) valutati
      mpz_t* As2; // A + s temporanei
      init_vector_mpz(&As2, intervals);

      init_vector_mpz(&evaluated_poly, intervals);

      max_fact += base_dim; // k + n

      pair * solutions = malloc(sizeof(pair) * base_dim);
      for(unsigned ii = 0; ii < base_dim; ++ii) {
	solutions[ii].sol1 = sol[ii].sol1;
	solutions[ii].sol2 = sol[ii].sol2;
      }
      
      if(startfrom != 0)
	for(unsigned ii = 0; ii < base_dim; ++ii) {
	  //printf("%d, sol1=%d, ", factor_base[ii], solutions[ii].sol1);
	  
	  //while(solutions[ii].sol1 < startfrom)
	  //	solutions[ii].sol1 += factor_base[ii];

	  unsigned int old1 = solutions[ii].sol1;
	  unsigned int old2 = solutions[ii].sol2;

	  //unsigned int f = (startfrom - solutions[ii].sol1) / factor_base[ii] + 1;
	  //solutions[ii].sol1 = solutions[ii].sol1 + f * factor_base[ii];
	  //f = (startfrom - solutions[ii].sol2) / factor_base[ii] + 1;
	  //solutions[ii].sol2 = solutions[ii].sol2 + f * factor_base[ii];


	  while(solutions[ii].sol1 < startfrom)
	    solutions[ii].sol1 += factor_base[ii];
	  while(solutions[ii].sol2 < startfrom)
	    solutions[ii].sol2 += factor_base[ii];
	    
	  //printf("%d, [old]sol1=%d, [new]sol1=%d [%d, %d]\n", factor_base[ii], old1, solutions[ii].sol1, startfrom, poly_val_num);
	}

      for(l = startfrom; (l <= poly_val_num - intervals) && go_on /*&& (flag == 0)*/; l += intervals) {
	
	//printf("%d\n", l);

	for(i = 0; i < ((intervals / N_BITS) + 1); ++i) {
	  set_matrix_l(is_used_expo2, 0, i, 0);
	}

	for(i = 0; i < intervals; ++i) {
	  mpz_add_ui(intermed, n_root, i + l);
	  mpz_mul(intermed, intermed, intermed);
	  mpz_sub(evaluated_poly[i], intermed, n);
			
	  mpz_add_ui(As2[i], n_root, i + l);
	}

	for(i = 0; i < base_dim && go_on; ++i) {
			
	  for(j = solutions[i].sol1; j < intervals + l && go_on; j += factor_base[i]) {

  //printf("j=%d\n", j);

	    while(mpz_divisible_ui_p(evaluated_poly[j - l], factor_base[i])) {

	      // Se non sono mai stati usati gli esponenti
	      if(get_k_i(is_used_expo2, 0, j - l) == 0) {
		for(k = 0; k < base_dim; ++k)
		  set_matrix(expo2, j - l, k, 0);
		set_k_i(is_used_expo2, 0, j - l, 1);
	      }

	      set_matrix(expo2, j - l, i, get_matrix(expo2, j - l, i) + 1); // ++exponents[j][i];
	      mpz_divexact_ui(evaluated_poly[j - l], evaluated_poly[j - l], factor_base[i]);
	      //printf("i=%d j=%d exp=%d\n", j - l, i, get_matrix(expo2, j - l, i));
	    }
	  }
	  solutions[i].sol1 = j; // Al prossimo giro ricominciamo da dove abbiamo finito

	  for(j = solutions[i].sol2; j < intervals + l && factor_base[i] != 2 && go_on; j += factor_base[i]) {

	    while(mpz_divisible_ui_p(evaluated_poly[j - l], factor_base[i])) {
	      // Se non sono mai stati usati gli esponenti
	      if(get_k_i(is_used_expo2, 0, j - l) == 0) {
		for(k = 0; k < base_dim; ++k)
		  set_matrix(expo2, j - l, k, 0);
		set_k_i(is_used_expo2, 0, j - l, 1);
	      }

	      set_matrix(expo2, j - l, i, get_matrix(expo2, j - l, i) + 1); // ++exponents[j][i];
	      mpz_divexact_ui(evaluated_poly[j - l], evaluated_poly[j - l], factor_base[i]);
	    }
	  }
	  solutions[i].sol2 = j; // Al prossimo giro ricominciamo da dove abbiamo finito
	}

	for(i = 0; i < intervals; ++i) {
	  if(mpz_cmp_ui(evaluated_poly[i], 1) == 0) {
	    ++fact_count;
	    //printf("thread %d) poly[%d] interval=%d\n", thread_id, i, intervals);
	    for(j = 0; j < base_dim; ++j)
	      buffer[j] = get_matrix(expo2, i, j);
	    buffer[j] = 0;
	    //printf("%d", buffer[j]);
	    //}
	    //gmp_printf(" - %Zd", As2[i]);
	    //printf("\n");
	    #pragma omp critical
	    {
	      MPI_Send(buffer, base_dim, MPI_UNSIGNED, 0, ROW_TAG, MPI_COMM_WORLD);
	      how_many_bytes = (mpz_sizeinbase(As2[i], 2) + 7) / 8;
	      *buffer_as = 0;
	      mpz_export(buffer_as, NULL, 1, 1, 1, 0, As2[i]);
	      MPI_Send(buffer_as, how_many_bytes, MPI_UNSIGNED_CHAR, 0, AS_TAG, MPI_COMM_WORLD);
	      
	      //printf("#%d) found!\n", thread_id);
	      //for(j = 0; j < base_dim; ++j)
	      //printf("%d", get_matrix(expo2, i, j));
	      //gmp_printf(" - %Zd\n", As2[i]);
	    }
	  }
	}

	char stop_signal;
	MPI_Request request;
	MPI_Status status;
	
	//MPI_Irecv(&stop_signal, 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &request);
	//MPI_Test(&request, &flag, &status);
	
	//printf("flag = %d\n", flag);

      }
      //printf("#%d)\n", thread_id);	      
    }

    // spedisco una riga del tipo [X X X X ... X 1] per indicare la terminazione
    unsigned int end_buffer[base_dim+1];
    end_buffer[base_dim] = 1;
    MPI_Send(end_buffer, base_dim+1, MPI_UNSIGNED, 0, ROW_TAG, MPI_COMM_WORLD);
    return 0;
}