# Brendan Philbin
# Simple Ising Model test w/ minimum energy analysis

declare -i TRIAL=1
declare -i M=100
declare -i B=5
declare -i MC_SEED=45
declare -i FERRO=0

for N in 16 32
do
    for J_SEED in 2 3 4 5 6
    do
        for SPIN_SEED in 7 8 9 10 11
        do
#            echo "Trial #$TRIAL : min_energy_output$TRIAL.csv : N = $N, M = $M, B = $B, J_SEED = $J_SEED, MC_SEED = $MC_SEED, SPIN_SEED = $SPIN_SEED, FERRO = $FERRO" >> receipt.txt
#            echo "" >> receipt.txt
            ./min_energy $N $M $B $J_SEED $MC_SEED $SPIN_SEED $FERRO $TRIAL
            ((TRIAL++))
        done
    done
done
