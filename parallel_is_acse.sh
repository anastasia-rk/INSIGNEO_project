#!/bin/bash

#$ -P acsehpc
#$ -q acsehpc.q
#$ -pe smp 16
#$ -l rmem=4G
#$ -cwd 
#$ -m bea                            
#$ -M a.kadochnikova@sheffield.ac.uk
#$ -o report_mf_uniform.txt
#$ -j y

module load apps/matlab/2019b
matlab -nodesktop -nosplash -r main_is_test