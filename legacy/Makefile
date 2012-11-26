%.o: %.f
	gfortran -g -O3 -pedantic -Wall -std=f95 -fimplicit-none -c -o $@ $^

geodc: geodc.o geodesicc.o
	cc -o $@ $^ -lm

geodf: geodf.o geodesicf.o
	gfortran -g -O3 -o $@ $^

geodesicc.o: geodesicc.h
geodc.o: geodesicc.h
