#include "../include/quadratic_sieve.h"

/* Ritorna un codice di errore oppure 0 */
unsigned int master(unsigned int base_dim, unsigned int max_fact, 
		    unsigned int** exponents, mpz_t * As,
		    int comm_size, unsigned int * n_fatt) {
  unsigned int fact_count = 0;
  unsigned int* buffer_exp;
  MPI_Status status1;
  MPI_Status status2;
  int count;
  int source;
  unsigned char buffer_As[BUFFER_DIM];

  /* Contatore degli slave che hanno terminato */
  unsigned int n_finished = 0;
  int si = 0;

  init_vector(& buffer_exp, base_dim+1);
  while(fact_count < max_fact + base_dim) {
    MPI_Recv(buffer_exp, base_dim+1, MPI_UNSIGNED, 
	     MPI_ANY_SOURCE, ROW_TAG, 
	     MPI_COMM_WORLD, &status1);
    
    //MPI_Get_count(&status2, MPI_UNSIGNED_CHAR, &count);
    //printf("Ho ricevuto %d\n", count);
    source = status1.MPI_SOURCE;
     
    //printf("%d) ", fact_count+1);
    //for(unsigned int k = 0; k < base_dim+1; ++k) {
    //  printf("%d", buffer_exp[k]);
    //}
    //printf("\n");

   MPI_Get_count(&status2, MPI_UNSIGNED_CHAR, &count);
   if(buffer_exp[base_dim] == 1) {
     ++n_finished;
     //printf("n_f=%d, comm=%d\n", n_finished, comm_size);
     if(n_finished >=  comm_size - 1) {
       *n_fatt = fact_count;
       return EVERYONE_FINISHED;
     }
   } else {
     
     /*
       if(count == 0) {
       //++n_finished;
       //if(n_finished >=  comm_size - 1) {
       //	*n_fatt = fact_count;
       //	printf("f=%u c=%u fatt=%u\n", n_finished, comm_size-1, *n_fatt);
     
       //return EVERYONE_FINISHED;
       //     }
       }*/
     
    
     //printf("Mm f=%d)\n", fact_count);
     for(unsigned int i = 0; i < base_dim; ++i) {
       set_matrix(exponents, fact_count, i, buffer_exp[i]);
       //  printf("%u", buffer_exp[i]);
     }
    
     MPI_Recv(buffer_As, BUFFER_DIM, MPI_UNSIGNED_CHAR, source, 
	      AS_TAG, MPI_COMM_WORLD, &status2);

     MPI_Get_count(&status2, MPI_UNSIGNED_CHAR, &count);
     mpz_import(As[fact_count], count, 1, 1, 1, 0, buffer_As);
    
     //gmp_printf("- %Zd\n", fact_count, As[fact_count]);
     //printf("\n");
     ++fact_count;
   }
  }
  // Spedisco '1' agli slave per indicare la terminazione
  //char stop_signal = '1';
  //MPI_Bcast(&stop_signal, 1, MPI_CHAR, 0, MPI_COMM_WORLD);
  //printf("M) Sending stop_signal\n");
  *n_fatt = fact_count;
  return 0;
}

unsigned long quadratic_sieve(mpz_t N, 
			      unsigned int n, 
			      unsigned int poly_val_num,
			      unsigned int max_fact,
			      unsigned int interval,
			      mpz_t m) {
  double t1, t2;
  
  int rank;
  int comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, & comm_size);

  /* Controllo con test di pseudoprimalità di rabin */
  if(mpz_probab_prime_p(N, 25)) {
    return NUM_PRIMO;
  }

  /* Radice intera di N */
  mpz_t s;
  mpz_init(s);
  mpz_sqrt(s, N); 

  t1 = MPI_Wtime();
  /* Individuazione primi in [2, n] */
  unsigned int * primes = malloc(sizeof(unsigned int) * n);  
  eratosthenes_sieve(primes, n);
  /* Compattiamo i numeri primi in primes */
  unsigned j = 0;
  for(int i = 2; i < n; ++i)
    if(primes[i] == 1) {
      primes[j++] = i;
    }

  unsigned int n_all_primes = j;
    
  /* Fattorizzazione eseguita da tutti, gli slave ritornano
     IM_A_SLAVE mentre il main il fattore */
  unsigned int simple_factor = trivial_fact(N, primes, n_all_primes);
  if(simple_factor != 0) {
    mpz_set_ui(m, simple_factor);
    return rank == 0 ? OK : IM_A_SLAVE;
  }
  /* Calcolo base di fattori e soluzioni dell'eq x^2 = N mod p */
  pair * solutions = malloc(sizeof(pair) * n_all_primes);
  unsigned int * factor_base = primes;

  unsigned n_primes = base_fattori(N, s, factor_base, solutions,
				  primes, n_all_primes);
  t2 = MPI_Wtime();
  double t_base = t2 - t1; 
  printf("#dimensione base di fattori: %d\n", n_primes);

  /* Vettore degli esponenti in Z */
  unsigned int ** exponents;
  /* Vettore degli (Ai + s) */
  mpz_t * As;
  /* Parte di crivello: troviamo le k+n fattorizzazioni complete */
  unsigned int n_fatt;

  t1 = MPI_Wtime();
  if(rank == 0){
    /* Inizializzazioni vettori */
    init_matrix(& exponents, n_primes + max_fact, n_primes);
    init_vector_mpz(& As, n_primes + max_fact);

    /* Procedura master che riceve le fatt. complete */
    master(n_primes, max_fact, exponents, As, comm_size, & n_fatt);
    //MPI_Abort(MPI_COMM_WORLD, 0);
  } else {
    unsigned int dom_decomp = poly_val_num / (comm_size-1);
    unsigned int final_point = dom_decomp * rank;
    unsigned int starting_point = final_point - dom_decomp;

    //printf("s=%d f=%d\n", starting_point, final_point);

    n_fatt = smart_sieve(N, factor_base, n_primes, solutions, 
			 final_point, max_fact, 
			 interval, starting_point);
    // per gli slave l'algoritmo termina qui
    // MPI_Finalize();
    return IM_A_SLAVE;
  }
  t2 = MPI_Wtime();
  double t_sieve = t2 - t1;
  printf("#numero fattorizzazioni complete trovate: %d\n", n_fatt);
  
  //for(unsigned int i = 0; i < n_fatt; ++i) {
  //   for(unsigned int k = 0; k < n_primes; ++k)
  //    printf("%d", get_matrix(exponents, i, k));
  //   gmp_printf(" - %Zd\n", As[i]);
  // }
  //printf("\n");
 
  t1 = MPI_Wtime();
  /* Matrice di esponenti in Z_2 organizzata a blocchi di bit */ 
  word ** M;
  /* Numero di blocchi di bit da utilizzare */
  unsigned long n_blocchi = n_primes / N_BITS + 1;
  /* Inizializzazione egli esponenti mod 2 */
  init_matrix_l(& M, n_fatt, n_blocchi);
  for(int i = 0; i < n_fatt; ++i)
    for(int j = 0; j < n_primes; ++j) {
      unsigned int a = get_matrix(exponents, i, j);
      set_k_i(M, i, j, a);
    }

  /* Vettore con le info (bit piu' a dx e num bit a 1) su M */
  struct row_stats * wt = malloc(sizeof(struct row_stats) * n_fatt);
  for(int i = 0; i < n_fatt; ++i)
    get_wt_k(M, i, n_primes, & wt[i]);

  /* In gauss gli esponenti sommati possono andare in overflow,
     li converto dunque in mpz */
  mpz_t ** exponents_mpz;
  mpz_t temp;
  mpz_init_set_ui(temp, 2);

  unsigned int a; 
  init_matrix_mpz(& exponents_mpz, n_fatt, n_primes);
  for(unsigned i = 0; i < n_fatt; ++i)
    for(unsigned j = 0; j < n_primes; ++j) {
      a = get_matrix(exponents, i, j);
      mpz_set_ui(temp, a);
      set_matrix_mpz(exponents_mpz, i, j, temp);
    }
  
  /* Eliminazione gaussiana */
  gaussian_elimination(exponents_mpz, M, As, N, 
		       n_fatt, n_primes, n_blocchi, wt);
  t2 = MPI_Wtime();
  double t_gauss = t2 - t1;

  /* In m ritorno un fattore non banale di N */
  unsigned int n_fact_non_banali = factorization(N, factor_base, 
						 M, exponents_mpz, 
						 As, wt, n_fatt, 
						 n_primes, m);
  
  printf("#time_base time_sieve time_gauss time_totale\n");
  printf("%.6f ", t_base);
  printf("%.6f ", t_sieve);
  printf("%.6f ", t_gauss);
  printf("%.6f\n", t_base + t_gauss + t_sieve);
  
  if(n_fact_non_banali > 0) {
    /* Pulizia della memoria */
    return OK;
  }
  else {
    /* Pulizia della memoria */
    return SOLO_FATTORIZZAZIONI_BANALI;
  }
}
   
