#!/bin/bash

NP=16
dir="ptt_fdr.a"

# Python backup in nscc1083
# export PATH=/home-fn/users/nscc1083/software/python.latest/bin/:$PATH
# export LD_LIBRARY_PATH=/home-fn/users/nscc1083/software/python.latest/lib/python2.7/site-packages/:$LD_LIBRARY_PATH
# export PYTHONPATH=/home-fn/users/nscc1083/software/python.latest/lib/python2.7/:/home-fn/users/nscc1083/software/python.latest/lib/python2.7/site-packages/:$PYTHONPATH

#export PATH=/home-fn/users/nscc1082/software/software/Python-2.7.9/bin/:$PATH
module load software/partitionfinder2

time python /home-user/software/partitionFinder/partitionfinder-2.1.1/PartitionFinderProtein.py -p $NP -v $dir

