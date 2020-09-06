all: pop_fixed

pop_fixed: fixed.cc aux.h replica.h jij.h
	g++ fixed.cc -o pop_fixed -g -fopenmp

clean:
	rm -f pop_fixed

cleanResults:
	rm -f ../results/*.csv
