#!/bin/bash

PATH_SENDER="../nn_models/standard-sc05-min3-nomax/cache_sars2/"
userName="vsteb002"
hostName="pearl.cis.fiu.edu"
PATH_RECEIVER="/aul/homes/vsteb002/emomis/data/06-2021-non-redundant/nn_models/standard-sc05-min3-nomax/"

scp -r $PATH_SENDER $userName"@"$hostName":"$PATH_RECEIVER
