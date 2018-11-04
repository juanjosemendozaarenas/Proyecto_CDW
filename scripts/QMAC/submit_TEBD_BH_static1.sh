#!/bin/bash
# Submit a job to the qmac cluste

#qsub ./TEBD_BH_static1.sh $*
#qsub -l h_vmem=5G ./TEBD_BH_static1.sh $*
qsub -pe threads 16 ./TEBD_BH_static1.sh $* 
#qsub -q all.q@dell-0-1.local ./TEBD_BH_static1.sh $* 
#qsub -q all.q@dell-0-48.local -l h_vmem=60G ./TEBD_BH_static1.sh $*
#qsub -q all.q@dell-0-48.local -l h_vmem=5G ./TEBD_BH_static1.sh $*
qstat  
