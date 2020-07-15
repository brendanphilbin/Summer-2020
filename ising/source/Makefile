all: min_energy ising_simple several_seeds annealing ising_pop

min_energy: ising_min_energy.cc
	g++ ising_min_energy.cc -o min_energy

ising_simple: ising_simple.cc
	g++ ising_simple.cc -o ising_simple

several_seeds: ising_several_seeds.cc
	g++ ising_several_seeds.cc -o several_seeds

annealing: ising_annealing.cc
	g++ ising_annealing.cc -o annealing

ising_pop: ising_pop.cc
	g++ ising_pop.cc -o ising_pop

cleanAll: cleanSource cleanResults

cleanSource:
	rm -f min_energy ising_simple several_seeds annealing ising_pop

cleanResults:
	rm -f ../simple/ising_output*.csv
	rm -f ../min_energy/min_energy_output*.csv ../min_energy/results.txt
	rm -f ../several_seeds/several_seeds_output* ../several_seeds/results.txt
	rm -f ../annealing/annealing_output* ../annealing/results.txt ../annealing/histogram*
	rm -f ../pop/pop_output* ../pop/results* ../pop/histogram*
