# Brendan Philbin
# Simple Ising Model test w/ minimum energy analysis

declare -i TRIAL=1
declare -i M=100
declare -i B=5
declare -i MC_SEED=45

for N in 16 32
do
    for J_SEED in 2 3 4 5 6
    do
        for SPIN_SEED in 7 8 9 10 11
        do
            ./min_energy -n $N -m $M -b $B -j $J_SEED -c $MC_SEED -s $SPIN_SEED -t $TRIAL
            ((TRIAL++))
        done
    done
done
