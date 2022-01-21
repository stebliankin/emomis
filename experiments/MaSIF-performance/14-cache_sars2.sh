#!/bin/bash

source utils/set_env_hpc.sh

PYTHONPATH=$(pwd):$PYTHONPATH; export PYTHONPATH

mkdir -p logs

echo "[$(date)] Start caching training data " #>> logs/caching_training_log.txt

python3 $PROJECT_DIR/src/masif_modules/cache_training_data.py nn_models.sc05_sars2.custom_params #1>> logs/caching_training_log.txt 2>logs/caching_err.txt

echo "[$(date)]Done " #>> logs/caching_training_log.txt

