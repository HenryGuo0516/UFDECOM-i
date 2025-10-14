#!/bin/sh
###BSUB -gpu "num=1:mode=exclusive_process"
#BSUB -gpu "num=1"
#BSUB -n 1
#BSUB -q gpu
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -J run_gpu
nvidia-smi > out

./UFDE_G


