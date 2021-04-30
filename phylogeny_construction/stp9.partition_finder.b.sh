#!/bin/bash

APP_NAME=AMD_small
NP=16
RUN="RAW"

dir="ptt_fdr.b"

# Python backup in nscc1083
# export PATH=/home-fn/users/nscc1083/software/python.latest/bin/:$PATH
# export LD_LIBRARY_PATH=/home-fn/users/nscc1083/software/python.latest/lib/python2.7/site-packages/:$LD_LIBRARY_PATH
# export PYTHONPATH=/home-fn/users/nscc1083/software/python.latest/lib/python2.7/:/home-fn/users/nscc1083/software/python.latest/lib/python2.7/site-packages/:$PYTHONPATH

export PATH=/home-fn/users/nscc1082/software/software/Python-2.7.9/bin/:$PATH

time python /home-fn/users/nscc1082/software/partitionfinder-master/PartitionFinderProtein.py -p 16 -v $dir


