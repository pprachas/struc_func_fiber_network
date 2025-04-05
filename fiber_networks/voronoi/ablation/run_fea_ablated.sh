#!/bin/bash -l

#$ -P lejlab2       # Specify the SCC project name you want to use
#$ -l h_rt=24:00:00  # Specify the hard time limit for the job
#$ -N run_ab  # Give job a name
#$ -j y              # Merge the error and output streams into a single file
#$ -m ea
#$ -t 1-2

module load miniconda/23.1.0
conda activate fenics

python3 run_fea_ablated.py $SGE_TASK_ID