#!/bin/bash


N=$1
M=$2
K=$3
id_batch=$4
ngr_batch=$5
eta=$6
steps_each=$7
tol=$8
nthreads=$9

donefile=${10}
project=${11}



for id in $(seq $id_batch $((id_batch + ngr_batch - 1)))
do
  cond=true
  done_str="${N} ${M} ${id} ${steps_each}" 
  while read line; do
    if [ "$line" = "$done_str" ]; then
      cond=false
      break
    fi
  done <$donefile

  if $cond; then

    errfile="Error_CDA-d_FMS_N_"$N"_M_"$M"_id_"$id"_stepsdec_"$steps_each".txt"

    "./CDA_decimation_FMS.out" $N $M $K $id $eta $steps_each $tol $nthreads > /dev/null 2> $errfile

    mv $errfile $project/Results/
    mv "CDA_decimation_FMS_dyn_K_"$K"_N_"$N"_M_"$M"_"**"_stepsdec_"$steps_each"_seed_"$id"_"* $project/Results/
    mv "CDA_decimation_FMS_final_K_"$K"_N_"$N"_M_"$M"_"**"_stepsdec_"$steps_each"_seed_"$id"_"* $project/Results/

    echo $done_str >> $donefile

  fi
done

exit 0

