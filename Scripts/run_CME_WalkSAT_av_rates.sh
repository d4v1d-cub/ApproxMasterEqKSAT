#!/bin/bash

mkdir Results
mkdir Programs

cd Programs

export LC_NUMERIC="en_US.UTF-8"

N=$1
K=$2
q=$3
tl=$4
tol=$5
seed0=$6
nsamples=$7
nthr=$8

alphamin=$9
alphamax=${10}
alphastep=${11}


for a in $(seq $alphamin $alphastep $alphamax)
do
for seed in $(seq $seed0 $((seed0 + nsamples - 1)))
do

./CME_WalkSAT_av_rates.out $N  $((N * a / 100)) $K $seed $q $tl $tol $nthr > "../Results/Out_CME_WalkSAT_av_rates_q_"$q"_alpha_"$a"_seed_"$seed".txt" 2> "../Results/Error_CME_WalkSAT_av_rates_q_"$q"_alpha_"$a"_seed_"$seed".txt" &

echo "done q="$q"  alpha="$a"	seed="$seed

done
done


exit 0

