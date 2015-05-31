


main: main.c factor_base.o legendre.o
	gcc main.c factor_base.o legendre.o -lgmp -o main

factor_base.o: factor_base.c
	gcc -c factor_base.c -std=c99

legendre.o: legendre.c
	gcc -c legendre.c -std=c99 -lgmp 

clean:
	rm *.o