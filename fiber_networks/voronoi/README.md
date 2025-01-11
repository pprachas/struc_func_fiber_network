# voronoi Directory
This directory contains code for analysis of our random fiber networks based on Voronoi diagrams. The code here is used to generate plots in our publication, but the exact format (i.e. locations of texts and legends) might be slightly different due to modifications done in Adobe Illustrator.

## List of code files
|File Name|Description|
----------|------------
|``crit_strain_voronoi.py``|Runs bisection solver and outputs critical strain transition point|
|``fea_voronoi.py``|Runs finite element solver and save results|
|``generate_voronoi.py``|Creates the mesh and graph for all random fiber networks used in this work|
|``mesh_refinement.py``|Performs mesh refinement|
|``mesh_voronoi.py``|Generates mesh for our random fiber networks|
|``res_convergence.py``|Runs our convergence analysis on the regularization term|
|``run_crit_strain.py``|Scripts to create directories and files to perform parameter sweep and save critical strain transition|
|``run_fea_voronoi.py``|Scripts to create directories and files to perform parameter sweep and save FEA results|

Note that bash scripts to run jobs in bulk are also provided for convenience.

## List of Directories
|Directory|Description|
----------|------------
|``plots``|Directory for code generating results plots|

Additional information on code in on plots can be found in ``plot`` directory.


