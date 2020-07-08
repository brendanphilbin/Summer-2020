TRIAL=1
B=5
R=100
N=64
J=2
S=3
C=4

for M in 100 200 400
do
    ./annealing -n $N -m $M -r $R -c $C -s $S -j $J -b $B -t $TRIAL
    ((TRIAL++))
    ((C++))
    ((S++))
done

C=4
S=3

for M in 100 200 400
do
    ./annealing -n $N -m $M -r $R -c $C -s $S -j $J -b $B -t $TRIAL -a
    ((TRIAL++))
    ((C++))
    ((S++))
done
