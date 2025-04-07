#!/bin/bash


export LC_NUMERIC="en_US.UTF-8"


pop_size=10000
K=3
seed_r=1
tol=1e-3
nthr=10
eps_c=1e-6

export OMP_NUM_THREADS=$nthr



eta=0.6
tl=20

for alpha in 2.98
do

./CDA1av_FMS_lpln_pop_dyn.out $pop_size $alpha $K $seed_r $eta $tl $tol $nthr $eps_c > "Out_K_"$K"_eta_"$eta"_alpha_"$alpha"_tl_"$tl"_popsize_"$pop_size".txt" 2> "Error_K_"$K"_eta_"$eta"_alpha_"$alpha"_tl_"$tl"_popsize_"$pop_size".txt" &

echo "eta="$eta "  alpha="$alpha"  sent"

done



eta=0.5
tl=20

for alpha in 3.17 3.20 3.23 3.27
do

./CDA1av_FMS_lpln_pop_dyn.out $pop_size $alpha $K $seed_r $eta $tl $tol $nthr $eps_c > "Out_K_"$K"_eta_"$eta"_alpha_"$alpha"_tl_"$tl"_popsize_"$pop_size".txt" 2> "Error_K_"$K"_eta_"$eta"_alpha_"$alpha"_tl_"$tl"_popsize_"$pop_size".txt" &

echo "eta="$eta "  alpha="$alpha"  sent"

done


eta=0.4
tl=20

for alpha in 3.55 3.59 3.63 3.66
do

./CDA1av_FMS_lpln_pop_dyn.out $pop_size $alpha $K $seed_r $eta $tl $tol $nthr $eps_c > "Out_K_"$K"_eta_"$eta"_alpha_"$alpha"_tl_"$tl"_popsize_"$pop_size".txt" 2> "Error_K_"$K"_eta_"$eta"_alpha_"$alpha"_tl_"$tl"_popsize_"$pop_size".txt" &

echo "eta="$eta "  alpha="$alpha"  sent"

done


eta=0.3
tl=20

for alpha in 4.12 4.14 4.16 4.18
do

./CDA1av_FMS_lpln_pop_dyn.out $pop_size $alpha $K $seed_r $eta $tl $tol $nthr $eps_c > "Out_K_"$K"_eta_"$eta"_alpha_"$alpha"_tl_"$tl"_popsize_"$pop_size".txt" 2> "Error_K_"$K"_eta_"$eta"_alpha_"$alpha"_tl_"$tl"_popsize_"$pop_size".txt" &

echo "eta="$eta "  alpha="$alpha"  sent"

done



eta=0.2
tl=20

for alpha in 5.00 5.05 5.10 5.15 5.20 5.25 5.30
do

./CDA1av_FMS_lpln_pop_dyn.out $pop_size $alpha $K $seed_r $eta $tl $tol $nthr $eps_c > "Out_K_"$K"_eta_"$eta"_alpha_"$alpha"_tl_"$tl"_popsize_"$pop_size".txt" 2> "Error_K_"$K"_eta_"$eta"_alpha_"$alpha"_tl_"$tl"_popsize_"$pop_size".txt" &

echo "eta="$eta "  alpha="$alpha"  sent"

done


wait

exit 0

