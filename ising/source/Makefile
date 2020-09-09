all: pop_fixed pop_fluc

pop_fixed: fixed.cc aux.h replica.h jij.h
	g++ fixed.cc -o pop_fixed -g -fopenmp

pop_fluc: fluc.cc aux.h replica.h jij.h
	g++ fluc.cc -o pop_fluc -g -fopenmp

clean:
	rm -f pop_fixed

cleanResults:
	rm -f ../results/*.csv
