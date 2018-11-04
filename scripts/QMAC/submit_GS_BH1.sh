#!/bin/bash
# Submit a job to the qmac cluste

#qsub ./GS_BH1.sh $*
#qsub -l h_vmem=5G ./GS_BH1.sh $*
#qsub -pe threads 16 ./GS_BH1.sh $*
qsub -q all.q@dell-0-5.local ./GS_BH1.sh $* 
qstat  
