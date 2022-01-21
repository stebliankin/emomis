#!/bin/bash
source utils/set_env_hpc.sh

N_RES=''

LOGS_DIR_CURR=${LOGS_DIR}/blast_${N_RES}res
mkdir -p $LOGS_DIR_CURR

curr_dir=$(pwd)
cd fasta/sabdab_antibodies/
for p in *; do
    cd $curr_dir
    echo $p
#    ./run_blast_one.sh $p
    sbatch -J blast \
                    --account=$SLURM_ACC \
                    --qos=$SLURM_QOS \
                    -p $SLURM_NODE_TYPE \
                    -N 1 \
                    -n 1 \
                    -o $LOGS_DIR_CURR"/stdout-${p}.txt" \
                    -e $LOGS_DIR_CURR"/stderr-${p}.txt" \
                    ./run_blast_one.sh $p fasta/masif_train_fasta/ blast_hits/ fasta/sabdab_antibodies/

done