TRIAL=1
N=32
R=50
J=2
S=3
C=4
B=5
for M in 1000 2000 4000
do
    ./annealing -n $N -m $M -r $R -j $J -s $S -c $C -b $B -t $TRIAL
    ((TRIAL++))
    ((C++))
    ((S++))
done

S=3
C=4

for M in 1000 2000 4000
do
    ./annealing -n $N -m $M -r $R -j $J -s $S -c $C -b $B -t $TRIAL -a
    ((TRIAL++))
    ((C++))
    ((S++))
done
