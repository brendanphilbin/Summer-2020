# Brendan Philbin
# Threaded fixed-population annealing test script for Ising Model

# N=512
R=36
M=10
K=10
B=5
J=3
C=2
S=4
L=5

# trial tracker
T=43

for N in 512 1024 2048 
do
    for G in 1 1 1 4 4 4
    do
        ./pop_fixed -n $N -r $R -m $M -k $K -b $B -t $T -j $J -c $C -s $S -l $L -g $G
        ((T=T+1))
    done
done
