#include <omp.h>
#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <stdio.h>
#include <unistd.h> /*contiene il prototipo di getopt() */
#include <string.h>


#include "../include/quadratic_sieve.h"


/*
Questo main possiede le opzioni
-n = numero da fattorizzare
-b = fattori primi nella base di ereatostene
-k = numero di fattorizzazioni addizionali da trovare
-i = grandezza dei blocchi degli slave
-r = numero di ripetizioni incrementali dell'algotitmo.


*/

/*
Parametri per la getopt()
*/
extern char *optarg;
extern int optind;
extern int opterr;


void show_usage(){
printf(" Utilizzo:");
printf("    -n numero da fattorizzare \n");
printf("    -b fattori primi nella base di ereatostene \n");
printf("    -k numero di fattorizzazioni addizionali da trovare \n");
printf("    -i grandezza dei blocchi degli slave \n");
printf("    -r numero di ripetizioni incrementali dell'algotitmo. \n");
}

int main(int argc, char ** argv) {


static char *optstring = "n:b:k:i:r:p:h";


if(argc == 1 ){
show_usage();
return 1;
}

char* Nn;                  // numero da fattorizzare
unsigned int n;             // numero di fattori in eratostene -n
unsigned int poly_val_num;  // intervallo per la valutazione del polinomio -M
unsigned int max_fact;      // numero di fattorizzazioni addizionali -k
unsigned int interval;      // numero di intervalli in cui spezzare sieve -i
unsigned int increments;    // numero di passi incrementali dell'algoritmo -r



int indice;
while ((indice=getopt(argc,argv,optstring))!=-1)
switch(indice){
    case 'n' :
        Nn = optarg;
        break;
    case 'b' :
        n =  atoi(optarg);
        break;
    case 'k' :
        max_fact = atoi(optarg);
        break;
    case 'i' :
        interval = atoi(optarg);
        break;
    case 'r' :
        increments = atoi(optarg);
        break;
    case 'M' :
        poly_val_num = atoi(optarg);
        break;
    case 'h':
        show_usage();
        return 1;
    case '?':
        show_usage();
        return 1;
    default :
        show_usage();
}







double t1, t2;
int retvalue = 0;

int rank;
MPI_Init(&argc, &argv);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);

mpz_t N;
mpz_init(N);

mpz_set_str(N, Nn, 10);


//gmp_printf("N: %Zd = %Zd * %Zd \n", N, P1, P2);



mpz_t N1;
mpz_t N2;
mpz_init(N1);
mpz_init(N2);

t1 = omp_get_wtime();
int status = quadratic_sieve(N, n, poly_val_num, max_fact, interval, N1);
//printf("ho finito %d stat = %d\n", rank, status);
t2 = omp_get_wtime();
double time = t2 - t1;

if(status == OK) {
mpz_divexact(N2, N, N1);
gmp_printf("#Fattorizzazione trovata: %Zd * %Zd = %Zd\n\n",
N1, N2, N);

printf("#N time\n");
gmp_printf("%Zd %.6f\n", N, time);

retvalue = 0;
}
else if(status == SOLO_FATTORIZZAZIONI_BANALI) {
printf("#Nessuna fattorizzazione non banale trovata\n\n");

//printf("#N time\n");
//gmp_printf("%Zd %.6f\n", N, time);

retvalue = 1;
}
else if(status == IM_A_SLAVE) {
//printf("#Slave %d chiude\n", rank);
retvalue = 0;
}

if(rank == 0)
MPI_Abort(MPI_COMM_WORLD, 0);
MPI_Finalize();
exit(0);

}
