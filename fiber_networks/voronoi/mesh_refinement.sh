#!/bin/bash -l

#$ -P lejlab2       # Specify the SCC project name you want to use
#$ -l h_rt=12:00:00  # Specify the hard time limit for the job
#$ -N mesh_refinement # Give job a name
#$ -j y              # Merge the error and output streams into a single file
#$ -m ea

module load miniconda/23.1.0
conda activate fenics

python3 mesh_refinement.py 