# voronoi Directory
This directory contains code for analysis of our random fiber networks based on Voronoi diagrams. The code here is used to generate plots in our publication, but the exact format (i.e. locations of texts and legends) might be slightly different due to modifications done in Adobe Illustrator.

## List of code files
|File Name|Description|
----------|------------
|crit_striain_voronoi.py||
|fea_voronoi.py||
|generate_voronoi.py|Creates the mesh and graph for all random fiber networks used in this work|
|mesh_refinement.py||
|mesh_voronoi.py||
|res_convergence.py||
|run_crit_strain.py||
|run_fea_voronoi.py||



## List of Directories
|Directory|Description|
----------|------------
|``discretized_chain``|FEA analysis on discretized chain|
|``plots``|Code to generate plots in results section of paper|
|``random_chain``|FEA analysis on random chain|
|``sinusoidal_chain``|FEA analysis on sinusoidal chain|
|``triangular_chain``|FEA analysis on triangular chain|

For all directories containing FEA analysis the files are arranged as:
|File Name|Description|
----------|------------
|``mesh_refinement.py``|Performs mesh refinement|
|``run_fea_xx.py``|Scripts to create directories and files to perform parameter sweep and save FEA results|
|``xx.py``|Runs bisection solver and outputs critical strain transition point |
|``run_xx.py``|Scripts to create directories and files to perform parameter sweep and save critical strain transition|

Note that bash scripts to run jobs in bulk are also provided for convenience.

Additional information on code in on plots can be found in ``plot`` directory.



