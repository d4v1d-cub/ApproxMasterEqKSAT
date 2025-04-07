#! /bin/bash

K=3
nsamples=1000
hist=1
path="../Results/"

graph_origin="inside"

tl=20
save_every=50

for N in 50000
do
for q in 0.0
do
for alpha in 60 80 100 120
do

M=$((alpha * N / 100))

nohup /home/machadod/julia-1.10.5/bin/julia ./walksat_greedy.JL $graph_origin $N $M $q $K $nsamples $hist $tl $path $save_every > "out_walksat-greedy_N_"$N"_q_"$q"_alpha_"$alpha".txt" 2> "err_walksat-greedy_N_"$N"_q_"$q"_alpha_"$alpha".txt" &
echo "N="$N"  q="$q"  alpha="$alpha""

done
done
done