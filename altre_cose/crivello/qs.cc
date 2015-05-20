/*
  Questo programma fa parte di una libreria messa a disposizione dei
  lettori del libro
  Alessandro Languasco & Alessandro Zaccagnini,
  Introduzione alla Crittografia, Ulrico Hoepli Editore, Milano, 2004.

  Scomposizione in fattori primi di un intero "n"
  Crivello quadratico di Pomerance
  Si veda il Paragrafo 6.6.1 del libro citato.
*/

#include <stdio.h>
#include <iostream>
#include <cmath>
#include "aritmetica.cc"

/*
  Il numero dei numeri primi precalcolati a disposizione
*/
static const
unsigned primi_prec = 40;

/*
  Il numero massimo dei numeri primi nella base di fattori
*/
static const
unsigned dim_bf = 15;

/*
  Il numero massimo di congruenze attese
*/
static const
unsigned max_num_congr = 5000;

/*
  Il numero di valori del polinomio ausiliario esplorati
*/
static const
unsigned max_num_val_poli = 15000;

/*
  Se "verbose = true", durante l'esecuzione vengono stampate
  informazioni dettagliate, e vengono creati alcuni files di dati.
*/
static const
bool verbose = false;

/*
  "primi_prec" numeri primi piccoli precalcolati per comodità.
  In un'applicazione reale si dovrebbe usare il Crivello di Eratostene.
*/
static const
long unsigned
primi[]={2,  3,  5,  7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53,
	 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113,
	 127, 131, 137, 139, 149, 151, 157, 163, 167, 173};

/*
  Alcuni file in cui vengono scritti i risultati delle varie fasi del
  calcolo.
*/
FILE *f_fb, *f_vec, *f_factors, *f_sol;

/*
  Prima del calcolo si verifica la condizione necessaria di Fermat
  "a^{n-1} = 1 mod n" per vari valori di "a".
  Se "n" soddisfa tutte queste condizioni allora è un "probabile primo".
  Piú precisamente, "n" è uno pseudoprimo in tutte le basi prime
  memorizzate in "primi[]"; quindi è probabilmente opportuno che sia
  sottoposto ad un criterio di primalità prima di tentare una sua
  scomposizione in fattori.
  NB: Le verifiche nelle basi composte sono superflue (se "n" è
  pseudoprimo in base 2 allora è pseudoprimo in base 4; se "n" è
  pseudoprimo in base 2 e 3, allora è pseudoprimo anche in base 6,
  ecc.) Viceversa, se "n" _non_ è pseudoprimo in base 2 o in base 3,
  il calcolo si interrompe senza raggiungere il valore j = 4, ecc.

  La funzione restituisce
  - "false" se trova un valore "a" nell'elenco dato sopra per cui
     a^(n-1) != 1 mod n,
  - "true" se a^(n-1) = 1 _per tutti_ gli "a" di cui sopra, cioè se
    "n" è uno pseudoprimo in tutte le basi esaminate.
*/
bool
criterio_di_fermat(long int n) {

  for (unsigned j = 0; j < primi_prec; ++j) {
    long unsigned p = primi[j];
    long i = potenze_mod_n(p, n - 1, n);
    // Se i != 1 allora "n" non è primo
    if (i != 1) {
      std::cout << "N non supera il criterio di Fermat in base " << p;
      std::cout << ": " << p << "^" << n - 1 << " = " << i;
      std::cout << " mod "<< n << std::endl;
      return(false);
    }
  }
  // "n" è un probabile primo
  return(true);
}

/*
  Determinazione della base di fattori, cioè dei primi "p" "piccoli"
  (quelli memorizzati nell'array "primi[]") per cui il simbolo di
  Legendre (n | p) = 1, cioè per cui l'equazione "x^2 = n mod p" ha
  soluzione.
  La funzione restituisce il numero di elementi della base di fattori.
  Le soluzioni dell'equazione "Q(A) = n mod p" (quelle che determinano
  le classi di congruenza da esplorare) sono memorizzate in
  "soluzioni[][]". I primi nella base di fattori sono copiati
  nell'array "primi_bf[]".
*/
long unsigned
base_di_fattori(long unsigned n, long unsigned s,
		long unsigned primi_bf[primi_prec],
		long unsigned soluzioni[primi_prec][2]) {

  // Il primo 2 appartiene sempre alla base di fattori
  primi_bf[0] = 2;
  // "k" contiene il numero di elementi della base di fattori trovati
  // finora
  long unsigned k = 1;
  for (unsigned i = 1; i < primi_prec; ++i) {
    long unsigned p = primi[i];
    long unsigned m = n % p;
    // Data la simmetria delle soluzioni di "x^2 = n mod p" quando "p"
    // è dispari, è sufficiente esaminare i valori di "x"
    // nell'intervallo [0, (p - 1) / 2]
    for (unsigned j = 0; 2 * j < p ; ++j)
      if ((j*j) % p == m) {
	// Se si entra qui vuol dire che "n" è un residuo quadratico
	// modulo "p", e quindi "p" fa parte della "base di fattori".
	// Le due soluzioni di "Q(A) = n mod p" si possono calcolare a
	// partire dalle soluzioni di "x^2 = n mod p".  Bisogna
	// accertarsi che la soluzione calcolata cada nell'intervallo
	// [0, p - 1], e per questo sono necessarie le operazioni di
	// modulo.
	soluzioni[k][0] = (p - ((s - j) % p)) % p;
	soluzioni[k][1] = (p - ((j + s) % p)) % p;
	// Memorizza il k-esimo elemento della base di fattori e
	// incrementa "k"
	primi_bf[k] = p;
	++k;
      }
  }
  if (verbose) {
    for (unsigned i = 0; i < k; ++i)
      printf("%5ld", primi_bf[i]);
    std::cout << std::endl;
  }

  // C'è una sola classe da esplorare modulo 2 durante la fase di
  // crivello, perché, se "N" è dispari, l'equazione "Q(A) = 0 mod 2"
  // diventa "(A + s) = 1 mod 2", la cui soluzione è "A = (s + 1) mod 2".
  soluzioni[0][0] = (s + 1) % 2;
  soluzioni[0][1] = (s + 1) % 2;
  if (verbose) {
    fprintf(f_sol, "&       &");
    for (unsigned i = 0; i < dim_bf; ++i)
      fprintf(f_sol, "&  %3ld  &", primi_bf[i]);
    fprintf(f_sol, "\\cr\n");
    for (int j = 0; j < 2; ++j) {
      fprintf(f_sol, "&       &");
      for (unsigned i = 0; i < dim_bf; ++i)
	fprintf(f_sol, "&  %3ld  &", soluzioni[i][j]);
      fprintf(f_sol, "\\cr\n");
    }
  }
  return(k);
}

/*
  Stampa le fattorizzazioni in varie forme in vari files ausiliari
*/
void
stampa_fattori(long int a, long int b, unsigned k,
	       long int c, long unsigned primi_bf[primi_prec],
	       long unsigned vettore[max_num_congr][dim_bf]) {

  fprintf(f_fb,      "A=%6ld; Q(A)=%11ld =", a, b);
  fprintf(f_factors, "& %6ld && %11ld &&", a, b);
/*    fprintf(f_vec,     "%5ld ", c + 1); */
  fprintf(f_vec,     "&%6ld &&%11ld && (", a, b);
  bool flag = false;
  for (unsigned i = 0; i < k; ++i) {
    if (vettore[c][i] > 0) {
      if (flag) {
	fprintf(f_fb, " *");
	fprintf(f_factors, "\\cdot");
      }
      else
	flag = true;
      fprintf(f_fb, "%3ld", primi_bf[i]);
      fprintf(f_factors, "%3ld", primi_bf[i]);
      if(vettore[c][i] > 1) {
	fprintf(f_fb,      "^%1ld", vettore[c][i]);
	fprintf(f_factors, "^%1ld", vettore[c][i]);
      }
      else {
	fprintf(f_fb, "  ");
	fprintf(f_factors, "  ");
      }
    }
    fprintf(f_vec,"%1ld", vettore[c][i]%2);
    if (i < k - 1)
      fprintf(f_vec, ",");
  }
  fprintf(f_fb, "\n");
  fprintf(f_factors, " & \\cr\n");
  // fprintf(f_vec," -- %5ld %3ld %3ld", c, wt[c][0], wt[c][1]);
  fprintf(f_vec,") & \\cr\n");
}

void
somma_mod_2(long int a, long int b, long int c,
	    long unsigned vettore_mod_2[max_num_congr][dim_bf]) {

  for (long i = 0; i < c; ++i) {
    long j = vettore_mod_2[a][i] + vettore_mod_2[b][i];
    vettore_mod_2[b][i] = j % 2;
  }
}

void
somma(long unsigned a, long unsigned b, long unsigned c,
      long unsigned vettore[max_num_congr][dim_bf]) {

  for (unsigned i = 0; i < c; ++i)
    vettore[b][i] = vettore[a][i] + vettore[b][i];
}

void
somma_bis(long unsigned a, long unsigned b, long unsigned c,
	  unsigned dep_bis[max_num_congr][200]) {

  for (unsigned i = 0; i < c; ++i)
    dep_bis[b][i] = (dep_bis[a][i] + dep_bis[b][i]) % 2;
}

/*
  Funzione ausiliaria
*/
void
get_wt(long unsigned a, long unsigned b,
       long unsigned vettore_mod_2[max_num_congr][dim_bf],
       long unsigned wt[max_num_congr][2]) {

  wt[a][0] = 0;
  wt[a][1] = b;
  for (long i = b - 1; i > -1; --i) {
    wt[a][0] += vettore_mod_2[a][i];
    if (vettore_mod_2[a][i])
      wt[a][1] = i;
  }
}

/*
  Questa funzione realizza il crivello.
  Sia "Q(A) = (A + s)^2 - N", dove "N" è il numero da scomporre in
  fattori, ed "s" è la sua radice quadrata intera (approssimata per
  difetto).
  "num_primi" è il numero di elementi nella base di fattori.
  1. Si memorizzano i valori di "Q(A)" per "A = 0", ...,
     "max_num_val_poli - 1" nel vettore "valori_pol_aus[]".
  2. Per ogni primo "p" nella "base di fattori" determinata
     dall'apposita funzione, si esplorano le 2 classi di congruenza
     mod p (una sola se "p = 2"), memorizzate nel vettore
     "soluzioni[][]" per cui "Q(A) = 0 mod p" se e solo se "A" è in
     una delle due classi, e si divide ripetutamente "Q(A)" per "p",
     finché ciò è possibile, memorizzando l'esponente di "p"
     corrispondente nel vettore "fattor_pol_aus[][]".
  3. I valori di "A" per cui "Q(A)" si fattorizza completamente sulla
     base di fattori, e le corrispondenti fattorizzazioni di "Q(A)",
     sono memorizzati per essere usati nella fase della ricerca delle
     dipendenze lineari.
  La funzione ritorna il numero di fattorizzazioni complete determinate.
*/
long unsigned
crivello(long unsigned n, long unsigned s,
	 long unsigned num_primi,
	 long unsigned primi_bf[primi_prec],
	 long unsigned soluzioni[primi_prec][2],
	 long unsigned vettore_mod_2[max_num_congr][dim_bf],
	 long unsigned wt[max_num_congr][2],
	 long unsigned vettore[max_num_congr][dim_bf],
	 long unsigned valori[max_num_congr],
	 long unsigned dep[max_num_congr][2],
	 unsigned dep_bis[max_num_congr][200]) {

  // Vettore in cui memorizziamo i valori del polinomio ausiliario
  long valori_pol_aus[max_num_val_poli] = {0};

  // Calcolo dei valori di "Q(i)"
  long unsigned i;
  for (i = 0; i < max_num_val_poli; ++i) {
    valori_pol_aus[i] = s * s - n + i * (i + 2 * s);
  }

  // Vettore in cui memorizziamo le fattorizzazioni degli interi nel
  // vettore "valori_pol_aus[]" rispetto alla base di fattori.
  unsigned fattor_pol_aus[max_num_val_poli][dim_bf] = {0};

  for (unsigned j = 0; j < num_primi; ++j) {
    unsigned p = primi_bf[j];
    // Vi sono 2 classi di congruenza da esplorare per tutti i primi
    // dispari, ma una sola per "p = 2"
    unsigned num_classi = 2;
    if (p == 2)
      num_classi = 1;
    for (unsigned k = 0; k < num_classi; ++k) {
      for (i = soluzioni[j][k]; i < max_num_val_poli; i += p) {
	long unsigned n = valori_pol_aus[i];
	while (n % p == 0) {
	  ++fattor_pol_aus[i][j];
	  n /= p;
	}
	valori_pol_aus[i] = n;
      }
    }
  }
  // Qui termina il crivello

  if (verbose) {
    std::cout << "Per questi valori il polinomio si fattorizza";
    std::cout << " completamente sulla base di fattori";
  }
  unsigned num_fatt_compl = 0;
  for (i = 0; i < max_num_val_poli; ++i)
    if (valori_pol_aus[i] == 1) {
      if (verbose)
	std::cout << i << "  ";
      // L'obiettivo di questa parte del codice è costruire la
      // matrice che contiene gli esponenti (rispetto ai primi nella
      // base di fattori) delle fattorizzazioni dei valori di "Q(A)"
      // che si fattorizzano completamente sulla base di fattori, e
      // un'altra matrice che contiene gli stessi esponenti ridotti
      // modulo 2. È su quest'ultima matrice che verrà eseguita
      // l'eliminazione di Gauss.

      // Per l'eliminazione di Gauss sono utili due ulteriori
      // quantità per ciascuna riga della matrice costruita come
      // detto sopra:
      // - il suo "peso" in bit, cioè il numero di bit = 1,
      //   memorizzato in "wt[][0]";
      // - la posizione del bit piú a sinistra in ciascuna riga,
      //   memorizzato in "wt[][1]".
      // Per questo motivo, ciascuna riga è scorsa da destra verso
      // sinistra e ad ogni bit = 1 viene incrementato il valore di
      // "wt[][0]" e aggiornato il valore di "wt[][1]".
      wt[num_fatt_compl][1] = num_primi;
      for (int j = num_primi - 1; j > -1; --j) {
	long unsigned l             = fattor_pol_aus[i][j] % 2;
	vettore      [num_fatt_compl][j] = fattor_pol_aus[i][j];
	vettore_mod_2[num_fatt_compl][j] = fattor_pol_aus[i][j] % 2;
	wt[num_fatt_compl][0] += l;
	if (l > 0)
	  wt[num_fatt_compl][1] = j;
	//    printf("%6ld %6ld\n", i, num_fatt_compl);
      }

      // Queste quantità serviranno per costruire la congruenza
      // "X^2 = Y^2 mod n" alla fine del crivello quadratico
      valori[num_fatt_compl] = i;
      dep[num_fatt_compl][0] = i + s;
      long unsigned q = s * s - n + i * (i + 2 * s);
      dep[num_fatt_compl][1] = q;
      dep_bis[num_fatt_compl][num_fatt_compl] = 1;
      //   printf("%8ld  %8ld\n", i, num_fatt_compl);
      if (verbose)
	stampa_fattori(i, q, num_primi, num_fatt_compl, primi_bf, vettore);
      ++num_fatt_compl;
    }
  if (verbose)
    std::cout << std::endl;
  std::cout << "Ho trovato " << num_fatt_compl << " fattorizzazioni complete";
  std::cout << std::endl << std::endl;
  return(num_fatt_compl);
}

/*
  Questa funzione realizza il crivello quadratico. Si noti che, nel
  caso sia in grado di produrre un fattore di "N", non c'è garanzia
  che questo fattore sia un numero primo, e quindi è opportuno
  chiamare ricorsivamente la stessa funzione sui fattori primi, fino a
  determinarne la fattorizzazione completa.
*/
bool
crivello_quadratico(long unsigned N) {

  std::cout << "Numero da fattorizzare: N = " << N << std::endl;
  // 1. Ricerca di eventuali fattori primi "piccoli" di "N"
  std::cout << "Passo 1: ricerca di eventuali fattori primi piccoli";
  std::cout << std::endl;
  for (unsigned i = 0; i < primi_prec; ++i) {
    long unsigned p = primi[i];
    if (N % p == 0) {
      unsigned alpha = 0;
      while (N % p == 0) {
	N /= p;
	++alpha;
      }
      std::cout << "Ho trovato il fattore " << p;
      if (alpha > 1)
	std::cout << "^" << alpha;
      std::cout << std::endl;
    }
  }
  long unsigned p = primi[primi_prec - 1];
  // Se abbiamo già rimosso tanti fattori primi, può darsi che la
  // parte non fattorizzata di "N" sia un numero primo, ed allora non
  // è necessario eseguire il resto del programma.
  if (p * p > N) {
    if (N > 1)
      std::cout << "Ho trovato il fattore " << N << std::endl;
    std::cout << "Fattorizzazione completata" << std::endl;
    return(true);
  }

  std::cout << "Numero da fattorizzare: N = " << N << std::endl;
  long s = radice_quadrata(N);
  std::cout << "Radice quadrata:        s = " << s << std::endl << std::endl;

  // 2. Verifica della condizione necessaria per la primalità
  std::cout << "Passo 2: verifica della condizione di Fermat" << std::endl;
  bool fermat = criterio_di_fermat(N);
  if (fermat)
    std::cout << N << " è un probabile primo" << std::endl;
  std::cout << std::endl;

  // 3. Determinazione dei primi piccoli "p" per cui "( N | p ) = 1",
  // cioè dei primi "p" per i quali l'equazione "x^2 = N mod p" ha
  // soluzione
  std::cout << "Passo 3: determinazione della base di fattori" << std::endl;

  // "primi_prec" coppie di eventuali soluzioni di "x^2 = N mod p"
  long unsigned soluzioni[primi_prec][2] = {0};

  // I numeri primi nella "base di fattori".
  long unsigned primi_bf[primi_prec] = {0};

  long unsigned num_primi = base_di_fattori(N, s, primi_bf, soluzioni);
  std::cout << "La base di fattori contiene " << num_primi << " numeri primi";

  // Permettiamo al massimo "dim_bf" numeri primi nella base di
  // fattori per due motivi: occupazione di memoria, e numero di
  // dipendenze lineari necessarie
  if (num_primi > dim_bf) {
    num_primi = dim_bf;
    std::cout << ". Consideriamo solo i primi " << dim_bf;
  }
  std::cout << std::endl << std::endl;

  // Array in cui si memorizza la scomposizione di "Q(A)" in fattori
  // primi appartenenti alla base di fattori.
  long unsigned vettore[max_num_congr][dim_bf] = {0};

  // Valori di "A" per cui "Q(A)" si fattorizza completamente sulla
  // base di fattori.
  long unsigned valori[max_num_congr] = {0};

  // Array in cui si memorizza la scomposizione di "Q(A)" in fattori
  // primi appartenenti alla base di fattori, con esponenti ridotti
  // mod 2
  long unsigned vettore_mod_2[max_num_congr][dim_bf] = {0};

  // Qui sono contenute le dipendenze lineari trovate alla fine del
  // crivello quadratico.
  long unsigned dep[max_num_congr][2] = {0};

  unsigned dep_bis[max_num_congr][200] = {0};

  // Matrici ausiliarie in uso nella ricerca delle dipendenze lineari
  // - wt[i][0] contiene il peso in bit (il numero di bit di valore 1)
  //   dell'i-esima riga;
  // - w[i][1] la posizione del bit = 1 piú a sinistra nell'i-esima
  //   riga.
  long unsigned wt[max_num_congr][2] = {0};

  // 4. Crivello
  std::cout << "Passo 4: il crivello" << std::endl;
  long unsigned fb = crivello(N, s, num_primi, primi_bf, soluzioni,
				   vettore_mod_2, wt, vettore, valori, dep,
				   dep_bis);

  // A questo punto
  // - "dep[j][0]" contiene il (j+1) esimo valore di "A"
  //   per cui "Q(A)" si fattorizza sulla base di fattori;
  // - "dep[j][1]" contiene "Q(A)";
  // - "wt[j][0]" contiene il "peso" della fattorizzazione di "Q(A)",
  //   cioè il numero di esponenti dispari nella fattorizzazione;
  // - "wt[j][1]" contiene la posizione del bit piú a sinistra = 1.

  // 5. algebra lineare: eliminazione di Gauss
  std::cout << "Passo 5: eliminazione di Gauss" << std::endl;
  for (unsigned i = 0; i < num_primi; ++i) {
    unsigned j;
    // Cerchiamo il primo vettore che abbia il bit piú a sinistra
    // esattamente nella posizione i-esima
    for (j = 0; (j < fb) & (wt[j][1] != i); ++j);
    // Usiamo questo vettore per eliminare un eventuale bit = 1 nella
    // posizione i-esima da tutti i vettori successivi. Teniamo
    // traccia di queste eliminazioni (e quindi delle dipendenze
    // lineari che troviamo) aggiornando i valori degli array
    // "dep[][]"
    for (unsigned k = j + 1; k < fb; ++k) {
      if (vettore_mod_2[k][i]) {
	somma_mod_2(j, k, num_primi, vettore_mod_2);
	somma(j, k, num_primi, vettore);
	somma_bis(j, k, fb, dep_bis);
	dep[k][0] = moltiplicazione_mod_n(dep[j][0], dep[k][0], N);
	dep[k][1] = moltiplicazione_mod_n(dep[j][1], dep[k][1], N);
	get_wt(k, num_primi, vettore_mod_2, wt);
      }
    }
  }

  long unsigned successi = 0;
  // I vettori che sono stati "eliminati" hanno peso = 0
  for (unsigned i = 0; i < fb; i++)
    if (wt[i][0] == 0) {
      long j = 1;
      long l;
      long unsigned m;
      for (unsigned k = 0; k < num_primi; ++k) {
	l = vettore[i][k] / 2;
	m = potenze_mod_n(primi_bf[k], l, N);
	j = moltiplicazione_mod_n(j, m, N);
      }
      long X = dep[i][0];
      long Y = j;
      printf("** %12ld^2 = %12ld^2 mod N;", X, Y);
      // 5. calcolo del MCD
      m = mcd(X + Y, N);
      l = N / m;
      printf("  N = %12ld * %12ld\n", m, l);
      if ((m > 1) & (m < N))
	++successi;
    }

  std::cout << "Fattorizzazione di N = " << N;
  std::cout << ". Numero di successi: " << successi << std::endl;
  if (verbose) {
    fprintf(f_vec, "\n\n");
    for (unsigned i = 0; i < fb; ++i)
      if (wt[i][0] == 0) {
	printf("%4u --- ", i);
	fprintf(f_vec, "& ");
	for (unsigned k = 0; k < fb; ++k)
	  if (dep_bis[i][k] > 0) {
	    printf("%3u,", k);
	    fprintf(f_vec, "\\vec v(%2ld) + ", valori[k]);
	  }
	printf("\n");
	fprintf(f_vec, "\\equiv 0 \\bmod 2 \\\\\n");
      }
  }
  std::cout << "Fine dell'esecuzione" << std::endl;
  return (successi > 0);
}

int
main() {

  if (verbose) {
    f_fb      = fopen("FactorBase.dat", "w");
    f_factors = fopen("Fattorizzazioni.dat", "w");
    f_sol     = fopen("Soluzioni.dat", "w");
    f_vec     = fopen("Vettori.dat", "w");
  }

/*
  Alcuni interi interessanti con cui provare il codice: l'asterisco
  indica i valori del parametro "dim_bf" per cui il crivello
  quadratico ha successo (e cioè scompone "N" in due fattori non
  banali, ma non necessariamente primi), fermi restando i valori di
  tutti gli altri parametri.

  N            = fattorizzazione completa                     dim_bf
                                                           10   15   20
  65537        = primo
  654079       = 593 * 1103                                 *    *    *
  67108865     = 5 * 53 * 157 * 1613                        *    *    *
  220424623    = 337 * 593 * 1103                           *    *    *
  536813567    = 8191 * 65537                                         *
  1073741825   = 5^2 * 13 * 41 * 61 * 1321 = 2^{30} + 1     *    *    *
  4294967297LL = 641 * 6700417 = 2^{32} + 1                 *    *    *
  8616460799LL = 89681 * 96079                                   *    *
  9049465682LL = 2 * 28879 * 156679                         *    *    *

  408704709982LL = 2 * 399601 * 511391
  Per quest'ultimo intero è necessario porre dim_bf = 40,
  max_num_val_poli = 40000.
*/

  long unsigned test[] = {65537, 654079, 67108865, 220424623,
			       536813567, 1073741825, 4294967297LL,
			       8616460799LL, 9049465682LL,
			       408704709982LL};

  for (unsigned i = 0; i < sizeof(test) / sizeof(long); ++i) {
    long unsigned N = test[i];
    bool status = crivello_quadratico(N);
    if (!status) {
      std::cout << "Il crivello quadratico non ha trovato fattorizzazioni";
      std::cout << std::endl;
    }
    std::cout << "Premi un tasto per continuare...";
    char c;
    std::cin >> c;
    std::cout << std::endl;
  }

  if (verbose) {
    fclose(f_fb);
    fclose(f_factors);
    fclose(f_sol);
    fclose(f_vec);
  }
}
