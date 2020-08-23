all: ising threaded

ising: main.cc anneal.h aux.h cxxopts.hpp fixed_anneal.h jij.h replica.h
	g++ -g main.cc -o ising

threaded: parallel.cc anneal.h aux.h cxxopts.hpp fixed_anneal.h jij.h replica.h
	g++ -g parallel.cc -o threaded -fopenmp

clean:
	rm -f ising threaded

cleanResults:
	rm -f ../results/*.csv
