#!/bin/bash

source utils/set_env_hpc.sh


PDB_SABDAB='5IQ9_B.fasta'

PDB_MASIF='1A0G_1.fasta'
HITS_DIR='emboss_hits'

mkdir -p $HITS_DIR

python emboss_needle_one.py fasta/sabdab_antibodies/${PDB_SABDAB} \
        fasta/masif_single_fasta/${PDB_MASIF} \
        ${HITS_DIR} \
        95

