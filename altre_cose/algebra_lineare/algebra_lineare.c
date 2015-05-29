#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define NUM_BLOCCHI 156
#define K 6000
#define NUM_PRIMI 128
#define N_BITS 64
#define TYPE unsigned long

typedef TYPE word;

struct row_stats {
  // bit piu a destra
  long unsigned b_dx;
  // num di bit a 1
  long unsigned n_bit;
};

struct row_stats wt[K];

// Ritorna l'i-mo bit della k-ma riga
unsigned int get_k_i(unsigned long M[][NUM_BLOCCHI], unsigned long k, 
		     unsigned long i) {
  unsigned long I = i / N_BITS;
  unsigned long n_shift = N_BITS - ((i % N_BITS ) + 1);

  //printf("i=%lu, I=%lu, n_shift=%lu, get=%lu, ", i, I, n_shift, (M[k][I] >> n_shift) & 1);
  //print_bits((M[k][I] >> n_shift));
  //printf(" \n");
  
  return (M[k][I] >> n_shift) & 1;
}

// Sostituisco all riga k la riga k + j:
//   k = k + j
// Utilizzo lo XOR bit a bit (che corrisponde)
// alla somma in modulo 2.
// Eseguo lo XOR tra ogni blocco dei vettori
unsigned long add_j_to_k(unsigned long M [][NUM_BLOCCHI], unsigned long k, 
			 unsigned long j, unsigned long n_blocchi) {
  for(unsigned long I = 0; I < n_blocchi; ++I)
    M[k][I] = M[k][I] ^ M[j][I];
}


void get_wt_k(unsigned long M[][NUM_BLOCCHI], unsigned long n_blocchi,
	      unsigned long k, struct row_stats * wt) {
  wt->b_dx = 128;
  wt->n_bit = 0;

  unsigned long i = 0;
  while(get_k_i(M, k, i) == 0 && i < (n_blocchi * N_BITS))
    ++i;

  if(i >= (n_blocchi * N_BITS))
    return;

  //printf("i=%lu get=%lu, ", i, get_k_i(M, k, i));

  wt->b_dx = i;

  for(i; i < (n_blocchi * N_BITS); ++i)
    if(get_k_i(M, k, i))
      wt->n_bit++;

  //printf("b_dx=%lu n_bit=%lu\n",  wt->b_dx,  wt->n_bit);
}

/*
void get_wt_k(unsigned long M[][NUM_BLOCCHI], unsigned long n_blocchi,
	      unsigned long k, struct row_stats * wt) {
  unsigned long b; // appoggio per copie locali
  unsigned long last_bit_pos = 0; // posizione ultimo bit a 1 nel blocco
  unsigned long bit_1_count = 0; // conta dei bit a 1
  unsigned long last_block_pos = 0; // blocco che contiene l'ultimo bit

  for(unsigned long I = n_blocchi; I > 0; --I) {
    b = M[k][I-1];
    //printf("%lu, I=%lu\n", M[k][I-1], I-1);
    for(unsigned i = 0; i < N_BITS; ++i) {
      if(b == 0)
	break;
      bit_1_count += b & 1; // sommo se il bit i-mo è a 1
      b = b >> 1;
      //printf("b=%lu, ", b);
      //printf("bit=%lu\n", bit_1_count);
      last_bit_pos = i;
      last_block_pos = I-1;
      //printf("last_bit=%lu, last_block=%lu\n", last_bit_pos, last_block_pos);
    }
  }

  printf("n_b=%lu b_dx=%lu I=%lu\n", n_blocchi*64-1, last_bit_pos, last_block_pos);

  (*wt).b_dx = (n_blocchi*64 - 1) - (last_bit_pos + last_block_pos*N_BITS);
  (*wt).n_bit = bit_1_count;
}
*/

void bit_gaussian_elimination_mod_2(unsigned long M[][NUM_BLOCCHI],
				    unsigned long n_row,
				    unsigned long n_col,
				    unsigned long n_blocks,
				    struct row_stats wt[]) {
  for(unsigned long i = 0; i < n_col; ++i) {
    unsigned long j;
    for(j = 0; j < n_row && wt[j].b_dx != i; ++j)
      //printf("wt[%d].b_dx=%d ==? %d\n", j, wt[j].b_dx, i)
      ;// avanzo e basta

    //printf("j=%d\n", j);

    for(unsigned k = j + 1; k < n_row; ++k) {
      //printf("wt[%d].b_dx=%d ==? %d\n", j, wt[j].b_dx, i);
      //printf("get=%d\n", get_k_i(M, k, i));
      //getchar();
      if(get_k_i(M, k, i)) { // il bit v(k)(i) deve essere a 1
	add_j_to_k(M, k, j, n_blocks); // v(k) = v(k) + v(j)
	//printf("add: %lu = %lu + %lu\n", k, k, j);

	// sommare le righe della matrice degli esponenti in Z
	
	// moltiplicare i Q(A)
	
	// aggiorno info su wt: bit piu' a destra e n bit a 1
	get_wt_k(M, n_blocks, k, & wt[k]);
      }
    }
  }
}

void print_bits(unsigned long a) {
  unsigned int bits[N_BITS];

  for(unsigned int i = 0; i < N_BITS; ++i)
    bits[i] = (a >> i) & 1U;

  for(int i = 63; i >= 0; --i)
    printf("%d", bits[i]);
}

void print_all(unsigned long M[][NUM_BLOCCHI], int righe){
  for(int i = 0; i < righe; ++i) {
    for(int j = 0; j < NUM_BLOCCHI; ++j) {
      print_bits(M[i][j]);
      printf(" ");
    }
    printf("\n");
  }
}

int main() {
  unsigned long M[K][NUM_BLOCCHI];

  /*       N_BITS          N_BITS
    1) 000 ... 001 000 ... 001
    2) 000 ... 000 000 ... 010
   */

  double t1 = clock();
  for(int i = 0; i < K; ++i)
    get_wt_k(M, NUM_BLOCCHI, i, & wt[i]);
  double t2 = clock();

  double t_set_up = (float) (t2 - t1)/CLOCKS_PER_SEC;

  //get_wt_k(M, 2, 1, & wt[1]);

  //for(int i = 0; i < 6; ++i)
  //printf("wt[].b_dx=%lu, wt[].n_bit=%lu\n", wt[i].b_dx, wt[i].n_bit);

  //print_all(M, K);

  //printf("\n\n");
  
  
  double t3 = clock();
  bit_gaussian_elimination_mod_2(M, K, NUM_BLOCCHI*N_BITS, NUM_BLOCCHI, wt);
  double t4 = clock();

  double t_gauss = (float) (t4 - t3)/CLOCKS_PER_SEC;

  printf("#time_gauss time_set_up time_totale\n");
  printf("%f ", t_gauss);
  printf("%f ", t_set_up);
  printf("%f\n", t_gauss + t_set_up);

  //print_all(M, K);

  //for(int i=0; i<2; i++){
  //print_bits(M[i][0]);
  //printf(" ");
  //print_bits(M[i][1]);
  //printf("\n");
    //get_wt_k(i, & wt);
    //printf("wt[].b_dx=%lu, wt[].n_bit=%lu\n", wt.b_dx, wt.n_bit); 
  //}
  
  //add_k_to_j(0, 1);

  //printf("\n");

  //print_bits(M[0][0]);
  //printf(" ");
  //print_bits(M[0][1]);
  //printf("\n");
  
  //for(int i = 0; i < 64; ++i)
    //get_k_i(0, i);
    //printf("%d", get_k_i(0, i));
	   
  //printf("\n");
}
