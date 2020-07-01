all: min_energy ising_simple

min_energy: ising_min_energy.cc
	g++ ising_min_energy.cc -o min_energy

ising_simple: ising_simple.cc
	g++ ising_simple.cc -o ising_simple

cleanAll: cleanSource cleanResults

cleanSource:
	rm -f min_energy ising_simple

cleanResults:
	rm -f ../simple/ising_output*.csv
	rm -f ../min_energy/min_energy_output*.csv ../results.txt
