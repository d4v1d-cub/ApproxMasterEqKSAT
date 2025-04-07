#! /bin/bash

N=$1
alpha=$2
K=$3
nsamples=$4
hist=$5
seed_fms=$6
iters_save=$7
tl=$8
eta=$9
path=${10}
print_every=${11}

M=$((alpha * N / 100))

for idumgraph in $(seq 1 $nsamples)
do

filehist=$path"FMS_KSAT_K_"$K"_eta_"$eta"_N_"$N"_alpha_"$alpha"_ngraphs_"$nsamples"_nhist_"$hist"_tl_"$tl"_idgraph_"$idumgraph".txt"
./Graph_to_CNF_input.out $N $M $K $idumgraph </dev/null 2>/dev/null | ./wfacwsat -FMS -restart $hist -seed $seed_fms -trace $iters_save -cutoff $((N * tl)) -print_every $print_every -fhist $filehist -noise $eta 100 > $path"out_FMS_K_"$K"_N_"$N"_eta_"$eta"_alpha_"$alpha"_idumgraph_"$idumgraph"_tl_"$tl"_dtN_"$iters_save".txt" 2>/dev/null
echo "sent  N="$N"  M="$M"  eta="$eta" idumgraph="$idumgraph

done

