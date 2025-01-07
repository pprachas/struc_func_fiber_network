#!/bin/bash -l

#$ -P lejlab2       # Specify the SCC project name you want to use
#$ -l h_rt=12:00:00   # Specify the hard time limit for the job
#$ -N run_random_HR           # Give job a name
#$ -j y               # Merge the error and output streams into a single file
#$ -m ea

# #$ -t 7-8 # number of job arrays


module load miniconda/23.1.0
conda activate fenics

python3 simulate_random_chain_HR.py 30 4 1000 3