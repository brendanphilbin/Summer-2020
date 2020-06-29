# Brendan Philbin
# Simple Ising Model test


# Infinite range Ising spin glass (J = +/- 1)
./ising_simple 4 100 0.1 3 5 1 0 > ising_output1.csv
./ising_simple 4 100 1 3 6 2 0 > ising_output2.csv
./ising_simple 4 100 10 3 7 3 0 > ising_output3.csv
./ising_simple 16 100 0.1 4 8 4 0 > ising_output4.csv
./ising_simple 16 100 1 4 9 5 0 > ising_output5.csv
./ising_simple 16 100 10 4 10 6 0 > ising_output6.csv
./ising_simple 64 100 0.1 5 11 7 0 > ising_output7.csv
./ising_simple 64 100 1 5 12 8 0 > ising_output8.csv
./ising_simple 64 100 10 5 13 9 0 > ising_output9.csv

# Infinite range ferromagnetic Ising model (J = +1)
./ising_simple 8 100 0.1 6 14 10 1 > ising_output10.csv
./ising_simple 8 100 1 6 15 11 1 > ising_output11.csv
./ising_simple 8 100 10 6 16 12 1 > ising_output12.csv
./ising_simple 16 100 0.1 7 17 13 1 > ising_output13.csv
./ising_simple 16 100 1 7 18 14 1 > ising_output14.csv
/ising_simple 16 100 10 7 19 15 1 > ising_output15.csv
./ising_simple 32 100 0.1 8 20 16 1 > ising_output16.csv
./ising_simple 32 100 1 8 21 17 1 > ising_output17.csv
./ising_simple 32 100 10 8 22 18 1 > ising_output18.csv
