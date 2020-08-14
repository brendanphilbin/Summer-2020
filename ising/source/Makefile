all: ising threaded shared

ising: main.cc anneal.h aux.h cxxopts.hpp fixed_anneal.h jij.h replica.h
	g++ -g main.cc -o ising

threaded: parallel.cc anneal.h aux.h cxxopts.hpp fixed_anneal.h jij.h replica.h
	g++ -g parallel.cc -o threaded -fopenmp

shared: shared.cc anneal.h aux.h cxxopts.hpp fixed_anneal.h jij.h replica.h
	g++ -g shared.cc -o shared -fopenmp

clean:
	rm -f ising threaded shared

cleanResults:
	rm -f ../results/*.csv
