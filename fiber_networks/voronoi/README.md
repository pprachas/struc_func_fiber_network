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

All the parameters in ``fea_voronoi.py`` and ``crit_strain_voronoi.py`` are the default ones. Some networks need different parameteres to converge. Those networks and parameters are:

|Network|Parameters|
--------|-----------
|n = 300, seed 17, kappa = 1e-5|init_c0 = 1e-6|
|n = 300, seed 17, kappa = 1e-6|init_c0 = 1e-6|
|n = 300 seed 1, kappa = 1e-4|init_c0 = 1e-6|
|n = 200 seed 9, kappa = 1e-6|init_c0 = 1e-5|
|n = 500 seed 6, kappa = 1e-6|init_c0 = 1e-5|

Note that the set of parameters to obtain convergence is not unique, but there are the ones that are used to obtain results in for our work.

Note that bash scripts to run jobs in bulk are also provided for convenience.



## List of Directories
|Directory|Description|
----------|------------
|``plots``|Directory for code generating results plots|
|``ablation``|Directory for code generating results and plots for fiber ablation studies|

Additional information on code in on plots can be found in ``plot`` directory.


