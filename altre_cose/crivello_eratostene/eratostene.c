/* crivello di eratostene
 * g++ eratostene.cc -o eratostene
 * use: ./eratostene -n NUMERO
 *      ./eratostene -h // HELP
 */

#include <errno.h>       /* error definitions and routines */
#include <stdlib.h>      /* C standard library */
#include <unistd.h>      /* unix standard library */
#include <stdio.h>	 /* standard I/O library */
#include <math.h>
#include <time.h>        /* clock() */

int just_times = 1;

long long unsigned n = 0;

void usage(void) {
  printf("Crivello di eratostene: calcola i primi tra 1 e n\n");
  printf("Usage:\n");
  printf("  primi [-h] [-n stop] [-v] \n");
  printf("  -h	   print this help\n");
  printf("  -n numero    numero \n");
  printf("  -v           solo n e tempi in output\n");
  exit(1);
}

int cli_read(int argc, char * argv[]){
  int i;
  if (argc==1)  usage();
  while ( (i = getopt(argc, argv, "hm:n:v")) != -1) {
    switch (i) {
    case 'h': 
      usage(); 
      return 1;
    case 'n':
      n = strtoll(optarg, NULL, 10); // numero
      break;
    case 'v':
      just_times = 0; // just times off
      break;
    default: 
      usage();
    }
    return 0;
  }
}

int main(int argc, char * argv[]){
  if(cli_read(argc, argv))
    return -1;

  if(just_times == 0)
    printf("#%llu\n", n);

  int * sieve = new int[n];
  
  // sieve[i] = 0 // composto
  // sieve[i] = 1 // primi
  for(long long unsigned i = 0; i <= n; i++)
    sieve[i] = 1;

  clock_t c1, c2;

  c1 = clock();

  for(long long unsigned i = 2; i <= n; i++){
    if(i*i > n)
      break;
    if(sieve[i] == 1)
      for(long long unsigned j = i; i*j <= n; j++)
	sieve[i*j] = 0;
  }

  c2 = clock();

  printf("#N time\n");
  printf("%lld %f\n", n, (float) (c2 - c1)/CLOCKS_PER_SEC);
  
  if(just_times == 0) {
    printf("# %d ", just_times);
    for(long long unsigned i = 2; i < n+1; ++i)
      if(sieve[i] == 1)
	printf("%llu ", i);
    printf("\n");
  }
  
}
