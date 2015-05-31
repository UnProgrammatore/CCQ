


main: main.c factor_base.o
	gcc main.c factor_base.o -lgmp -o main

factor_base.o: factor_base.c legendre.o
	gcc -c factor_base.c legendre.o -std=c99

legendre.o: legendre.c
	gcc -c legendre.c -std=c99 -lgmp 

clean:
	rm *.o