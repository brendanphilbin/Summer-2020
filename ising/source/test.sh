# Brendan Philbin
# Threaded fixed-population annealing test script for Ising Model

N=32
R=100
M=10
K=10
B=10
J=2
C=3
S=4
L=5

for T in 1 2 3 4
do
    ./threaded -n $N -r $R -m $M -k $K -b $B -t $T -j $J -c $C -s $S -l $L -g $T
done
