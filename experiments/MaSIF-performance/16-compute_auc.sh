#!/bin/bash

# RUN INSIDE SINGULARITY

#source utils/set_env_hpc.sh

#PYTHONPATH=$(pwd):$PYTHONPATH; export PYTHONPATH

mkdir -p logs

#echo "[$(date)] Start caching training data " >> logs/caching_training_log.txt

#python3 compute_AUC.py nn_models.sc05.custom_params #1>> logs/caching_training_log.txt 2>logs/caching_err.txt
#python3 compute_AUC.py nn_models.standard-sc05-min3-nomax.custom_params #1>> logs/caching_training_log.txt 2>logs/caching_err.txt

#python3 compute_AUC.py nn_models.sc05.custom_params,nn_models.standard-sc05-min3-nomax.custom_params,nn_models.sars2_masif_original.custom_params,nn_models.sars2_cache.custom_params

python3 compute_AUC.py nn_models.sc05.custom_params,nn_models.sars2_masif_original.custom_params
#python3 compute_AUC.py nn_models.sc05.custom_params

#echo "[$(date)]Done " >> logs/caching_training_log.txt

