# Test script for population annealing

N=64
R=100
M=10
K=100
B=10
J=1
C=2
S=3
L=4

for T in 1 2 3 4 5
do
    ./ising -n $N -r $R -m $M -k $K -b $B -j $J -c $C -s $S -l $L -t $T
    ((C++))
    ((S++))
    ((L++))
done
