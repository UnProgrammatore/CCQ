/*
  Questo programma fa parte di una libreria messa a disposizione dei
  lettori del libro
  Alessandro Languasco & Alessandro Zaccagnini,
  Introduzione alla Crittografia, Ulrico Hoepli Editore, Milano, 2004.
*/

#include <cmath>

/*
  Calcolo del massimo comun divisore con l'algoritmo di Euclide.
  Languasco & Zaccagnini, Introduzione alla Crittografia, Hoepli.
  Paragrafo 6.2
*/
long
mcd(long a, long b) {

  // I valori assoluti assicurano che i calcoli contengono solo numeri
  // non negativi
  long unsigned m = std::abs(a);
  long unsigned n = std::abs(b);
  long unsigned r = n;
  while(r != 0) {
    r = m % n;
    m = n;
    n = r;
  }
  return(m);
}

/*
  Calcolo del massimo comun divisore con l'algoritmo di Euclide esteso
  ritorna due valori "lambda" e "mi" in modo che sia soddisfatta
  l'uguaglianza "mcd = lambda * a + mi * b".
  Languasco & Zaccagnini, Introduzione alla Crittografia, Hoepli.
  Paragrafo 6.2
*/

long
mcd_est(long a, long b, long &lambda, long &mi) {

  long a1 = 1;
  long a2 = 0;
  long b1 = 0;
  long b2 = 1;
  long m = a;
  long n = b;
  long r = n;
  while(r != 0) {
    long q = m / n;
    r = m % n;
    long a3 = a1 - q * a2;
    long b3 = b1 - q * b2;
    a1 = a2;
    a2 = a3;
    b1 = b2;
    b2 = b3;
    m = n;
    n = r;
  }
  lambda = a1;
  mi = b1;
  return(m);
}

/*
  Moltiplicazione "a * m mod n" usando il procedimento del raddoppiamento
  e dimezzamento dei due fattori.
  Languasco & Zaccagnini, Introduzione alla Crittografia, Hoepli.
  Paragrafo 6.9.1
*/

long
moltiplicazione_mod_n(long a, long m, long n) {

  long S = 0;
  long M = m;
  long A = a;
  long r, q;
  do {
    r = M % 2;
    q = M / 2;
    if (r == 1) {
      S += A;
      S  = S % n;
    }
    A += A;
    A  = A % n;
    M  = q;
  } while(q > 0);
  return(S);
}

/*
  Calcolo di "a^m mod n" con il metodo dei quadrati ripetuti.
  Languasco & Zaccagnini, Introduzione alla Crittografia, Hoepli.
  Paragrafo 6.9.2
*/

long
potenze_mod_n(long a, long m, long n) {

  long P = 1;
  long M = m;
  long A = a;
  long r;
  long q;
  do {
    r = M % 2;
    q = M / 2;
    if (r == 1)
      P = moltiplicazione_mod_n(P,A,n);
    A = moltiplicazione_mod_n(A,A,n);
    M  = q;
  } while(q > 0);
  return(P);
}

/*
  Calcolo del quoziente "q" e del resto "r" della divisione di "N" per
  "M" sfruttando la divisione per "m = M - d", dove "d" è piccolo in
  valore assoluto (supponendo che, per qualche motivo, dividere per
  "m" sia efficiente).
  Languasco & Zaccagnini, Introduzione alla Crittografia, Hoepli.
  Paragrafo 6.9.4
*/

void
divisione(long N, long M, long d, long& q,
	  long& r) {

  long m = M - d;
  long max = d > 0 ? M : m;
  q = N / M;
  r = d * q + (N % M);
  long qq;

  // Ripeti il ciclo finché il resto non è nell'intervallo corretto
  do {
    qq = r / M;
    if (r < 0)
      --qq;
    long rr = r - M * qq;
    q += qq;
    r = d * qq + rr;
  } while ((r < 0) || (r >= max));

  // Correggi il risultato, se necessario
  while (r >= m) {
    r -= m;
    ++q;
  }
}

/*
  Radice quadrata intera.
  Dato il numero intero non negativo "n", restituisce il massimo
  intero non negativo "m" tale che "m^2 <= n", cioè la radice
  quadrata intera approssimata per difetto.
*/
long
radice_quadrata(long n) {

  if (n < 0)
    return(-1);
  long m = int(sqrt((double) n));
  return(m);
}

