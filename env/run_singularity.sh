#!/bin/bash


#cuda
#singularity shell --bind /opt,/disk,/usr/local/cuda/,/usr/lib64/,/usr/bin/ masif_latest.sif
# module load singularity-3.5.3
#singularity shell --bind /usr/local/cuda/,/usr/lib64/,/usr/bin/ emomis.sif

#module load singularity-3.5.3
singularity  shell emomis.sif

