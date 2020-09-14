#!/bin/bash

#$ -P acsehpc
#$ -q acsehpc.q
#$ -l h_rt=96:00:00
#$ -pe smp 16
#$ -l rmem=4G
#$ -cwd 
#$ -m bea                            
#$ -M a.kadochnikova@sheffield.ac.uk
#$ -o report_is_uniform_21_40.txt
#$ -j y

module load apps/matlab/2019b
matlab -nodesktop -nosplash -r main_is_41_60
