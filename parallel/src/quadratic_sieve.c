#include "../include/quadratic_sieve.h"

unsigned int master(unsigned int base_dim, unsigned int max_fact, 
		    unsigned int** exponents, mpz_t * As,
		    int comm_size, unsigned int print_fact) {

  unsigned int fact_count = 0;

  MPI_Status status;

  int count;
  int source;
  
  /* Buffer per ricevere gli esponenti */
  unsigned int* buffer_exp;
  /* Buffer per ricevere (A + s) */
  unsigned char buffer_As[BUFFER_DIM];
  init_vector(& buffer_exp, base_dim);

  double t1 = MPI_Wtime();
  double t2;

  while(fact_count < max_fact + base_dim) {
    /* Ricevo il vettore di esponenti */
    MPI_Recv(buffer_exp, base_dim, MPI_UNSIGNED,
	     MPI_ANY_SOURCE, ROW_TAG, 
	     MPI_COMM_WORLD, &status);
    source = status.MPI_SOURCE;

    if(fact_count % print_fact == 0) {
      t2 = MPI_Wtime() - t1;
      
      printf("#M) %d/%d in %.6f seconds\n", fact_count, max_fact + base_dim, t2);
    }
    
    for(unsigned int i = 0; i < base_dim; ++i) 
      set_matrix(exponents, fact_count, i, buffer_exp[i]);
    
    /* Ricevo l'mpz contenente (A + s) */
    MPI_Recv(buffer_As, BUFFER_DIM, MPI_UNSIGNED_CHAR, source, 
	     AS_TAG, MPI_COMM_WORLD, &status);
    MPI_Get_count(&status, MPI_UNSIGNED_CHAR, &count);
    mpz_import(As[fact_count], count, 1, 1, 1, 0, buffer_As);
    
    ++fact_count;
  }
  
  /* Spedisco '1' agli slave per indicare la terminazione */
  char stop_signal = '1';
  for(unsigned int i = 1; i < comm_size; ++i)
    MPI_Send(&stop_signal, 1, MPI_CHAR, i, 0, MPI_COMM_WORLD);
  
  printf("#M) Sending stop_signal\n");
  return fact_count;
}

unsigned long quadratic_sieve(mpz_t N, 
			      unsigned int n, 
			      unsigned interval,
			      unsigned int max_fact,
			      unsigned int block_size,
			      mpz_t m,
			      unsigned int print_fact) {
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
  if(rank == 0)
    printf("#Dimensione base di fattori: %d\n", n_primes);

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
    n_fatt = master(n_primes, max_fact, exponents, As, comm_size, print_fact);
  } else {
    mpz_t begin;
    mpz_init(begin);
    mpz_t counter;
    mpz_init(counter);

    mpz_set_ui(begin, interval * (rank - 1));
  
    //gmp_printf("%d) begin=%Zd interval=%d\n", rank, begin, interval);

    int stop_flag = 0;
    do {
      //gmp_printf("\t%d) [%Zd, %Zd+%d] - (flag=%d)\n", rank, begin, begin, interval, flag);
      stop_flag = smart_sieve(N, factor_base, n_primes, solutions,
		  begin, interval,
		  block_size, max_fact);
      mpz_add_ui(begin, begin, interval * (comm_size-1));
    } while(!stop_flag);

    printf("#%d) Termina\n", rank);

    return IM_A_SLAVE;
  }
  t2 = MPI_Wtime();
  double t_sieve = t2 - t1;
  printf("#Numero fattorizzazioni complete trovate: %d\n", n_fatt);
 
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
    return OK;
  }
  else {
    return SOLO_FATTORIZZAZIONI_BANALI;
  }
}
   
