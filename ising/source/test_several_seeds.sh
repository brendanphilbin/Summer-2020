# Test script for several replicas

# ./several_seeds -n 16 -m 100 -r 5 -t 1 -j 9 -c 2,3,4,5,6 -s 6,7,8,9,10 -b 10
# ./several_seeds -n 16 -m 100 -r 5 -t 2 -j 10 -c 3,4,5,6,7 -s 7,8,9,10,11 -b 10

TRIAL=1
M=100
B=5
R=5
N=64
J=10

for C in 1 6 11
do
    ./several_seeds -n $N -m $M -b $B -r $R -j $J -c $(($C)),$(($C+1)),$(($C+2)),$(($C+3)),$(($C+4)) -s $(($C+5)),$(($C+6)),$(($C+7)),$(($C+8)),$(($C+9)) -t $TRIAL
    ((TRIAL++))
done
