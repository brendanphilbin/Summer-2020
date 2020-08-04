all: ising 

ising: main.cc anneal.h aux.h cxxopts.hpp fixed_anneal.h jij.h replica.h
	g++ -g main.cc -o ising

clean:
	rm -f ising

cleanResults:
	rm -f ../results/*.csv
