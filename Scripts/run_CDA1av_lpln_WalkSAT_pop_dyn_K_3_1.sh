#!/bin/bash


export LC_NUMERIC="en_US.UTF-8"


pop_size=10000
K=3
seed_r=1
tol=1e-3
nthr=10
eps_c=1e-6

export OMP_NUM_THREADS=$nthr


q=0.2
tl=20

for alpha in 2.64 2.65
do

./CDA1av_WalkSAT_lpln_pop_dyn.out $pop_size $alpha $K $seed_r $q $tl $tol $nthr $eps_c > "Out_CDA1av_WalkSAT_K_"$K"_q_"$q"_alpha_"$alpha"_tl_"$tl"_popsize_"$pop_size".txt" 2> "Error_CDA1av_WalkSAT_K_"$K"_q_"$q"_alpha_"$alpha"_tl_"$tl"_popsize_"$pop_size".txt" &

echo "q="$q "  alpha="$alpha"  sent"

done


q=0.1
tl=20

for alpha in 2.67 2.68 2.69 2.70
do

./CDA1av_WalkSAT_lpln_pop_dyn.out $pop_size $alpha $K $seed_r $q $tl $tol $nthr $eps_c > "Out_CDA1av_WalkSAT_K_"$K"_q_"$q"_alpha_"$alpha"_tl_"$tl"_popsize_"$pop_size".txt" 2> "Error_CDA1av_WalkSAT_K_"$K"_q_"$q"_alpha_"$alpha"_tl_"$tl"_popsize_"$pop_size".txt" &

echo "q="$q "  alpha="$alpha"  sent"

done


exit 0

