#!/bin/bash

EMoMiS_PATH=$(pwd)/../../

#singularity exec --home /lclhome/${USER} $EMoMiS_PATH/env/emomis-singularity_v.2.5.2.sif python $EMoMiS_PATH/src/EMoMiS.py

singularity exec $EMoMiS_PATH/env/emomis.sif python $EMoMiS_PATH/src/EMoMiS.py --config config_variants.py --skip_blast --skip_sabdab
