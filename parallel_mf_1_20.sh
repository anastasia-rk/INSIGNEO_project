#!/bin/bash

#$ -P acsehpc
#$ -q acsehpc.q
#$ -l h_rt=96:00:00
#$ -pe smp 12
#$ -l rmem=4G
#$ -cwd 
#$ -m bea                            
#$ -M a.kadochnikova@sheffield.ac.uk
#$ -o report_mf_uniform_1_20.txt
#$ -j y

module load apps/matlab/2019b
matlab -nodesktop -nosplash -r main_mf_1_20