#!/bin/bash

EMoMiS_PATH=$(pwd)/../../

#singularity exec --home /lclhome/${USER} $EMoMiS_PATH/env/emomis-singularity_v.2.5.2.sif python $EMoMiS_PATH/src/EMoMiS.py --reverse

singularity exec $EMoMiS_PATH/env/emomis.sif python $EMoMiS_PATH/src/EMoMiS.py
