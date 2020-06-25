# Brendan Philbin
# Simple Ising Model test

# All J values +1
./ising_simple_j1 4 100 0.1 3 5 > ising_output1.csv
./ising_simple_j1 4 100 1 3 6 > ising_output2.csv
./ising_simple_j1 4 100 10 3 7 > ising_output3.csv
./ising_simple_j1 16 100 0.1 4 8 > ising_output4.csv
./ising_simple_j1 16 100 1 4 9 > ising_output5.csv
./ising_simple_j1 16 100 10 4 10 > ising_output6.csv
./ising_simple_j1 64 100 0.1 5 11 > ising_output7.csv
./ising_simple_j1 64 100 1 5 12 > ising_output8.csv
./ising_simple_j1 64 100 10 5 13 > ising_output9.csv

# Random J values +/- 1
./ising_simple_jrandom 8 100 0.1 6 14 > ising_output10.csv
./ising_simple_jrandom 8 100 1 6 15 > ising_output11.csv
./ising_simple_jrandom 8 100 10 6 16 > ising_output12.csv
./ising_simple_jrandom 16 100 0.1 7 17 > ising_output13.csv
./ising_simple_jrandom 16 100 1 7 18 > ising_output14.csv
./ising_simple_jrandom 16 100 10 7 19 > ising_output15.csv
./ising_simple_jrandom 32 100 0.1 8 20 > ising_output16.csv
./ising_simple_jrandom 32 100 1 8 21 > ising_output17.csv
./ising_simple_jrandom 32 100 10 8 22 > ising_output18.csv
