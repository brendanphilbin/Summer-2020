# Brendan Philbin
# Threaded fixed-population annealing test script for Ising Model

N=128
R=12
M=10
K=10
B=10
J=2
C=3
S=4
L=5

for T in 1 2 3 4
do
    ./pop_fixed -n $N -r $R -m $M -k $K -b $B -t $(($T+4)) -j $J -c $C -s $S -l $L -g $T
done
