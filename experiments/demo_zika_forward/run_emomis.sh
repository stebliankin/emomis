#!/bin/bash

EMoMiS_PATH=$(pwd)/../../

singularity exec $EMoMiS_PATH/env/emomis.sif python $EMoMiS_PATH/src/EMoMiS.py --config config_zika.py
