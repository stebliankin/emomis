#!/bin/bash
source utils/set_env_hpc.sh

LOGS_DIR=${LOGS_DIR}/emboss

mkdir -p $LOGS_DIR

curr_dir=$(pwd)

BLAST_HIT_DIR='blast_hits/'

cd $BLAST_HIT_DIR

for p in *; do
    cd $curr_dir
    echo $p
#    ./run_blast_one.sh $p
    sbatch -J emboss \
                    --account=$SLURM_ACC \
                    --qos=$SLURM_QOS \
                    -p $SLURM_NODE_TYPE \
                    -N 1 \
                    -n 1 \
                    -o $LOGS_DIR"/stdout-${p}.txt" \
                    -e $LOGS_DIR"/stderr-${p}.txt" \
                    ./pairwise_alignment_one.sh $p
done