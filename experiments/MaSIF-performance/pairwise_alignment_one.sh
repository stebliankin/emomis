#!/bin/bash

source utils/set_env_hpc.sh


BLAST_HIT=$1
SABDAB_FASTA_DIR='fasta/sabdab_antibodies/'
MASIF_FASTA_DIR='fasta/masif_single_fasta/'

HITS_DIR='emboss_hits'
IDENTITY_THRESHOLD=95
mkdir -p $HITS_DIR

python emboss_needle_one.py $BLAST_HIT \
    ${SABDAB_FASTA_DIR} \
    $MASIF_FASTA_DIR \
    ${HITS_DIR} \
    $IDENTITY_THRESHOLD

done
