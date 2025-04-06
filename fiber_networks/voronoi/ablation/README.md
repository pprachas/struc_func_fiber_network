# ablation Directory
This directory contains code for analysis of our random fiber networks based on Voronoi diagrams. The code here is used to generate plots in our publication, but the exact format (i.e. locations of texts and legends) might be slightly different due to modifications done in Adobe Illustrator.

## List of code files
|File Name|Description|
----------|------------
|``abalted_graph.py``|Plots ablated graph|
|``check_simulation.py``|Checks for simulations that did not converge|
|``compare ablated results.py``|plots normalized force vs. applied strain|
|``compare_relative_path_length.py``|compares new shortest oath with original second hsortest path|
|``generate_ablated_mesh.py``|Runs fea simulations of ablated fiber network|

All the parameters in ``fea_voronoi.py`` and ``crit_strain_voronoi.py`` are the default ones. Some networks need different parameteres to converge. Those networks and parameters are:

|Network|Parameters|
--------|-----------
|n=300 seed 1 edge 9| init_c0 = 1e-6|
|n=500 seed 1 edge 7| init_c0 = 1e-6|
|n=500 seed 1 edge 9| init_c0 = 1e-5|
|n=500 seed 1 edge 16| init_c0 = 1e-6|
|n=500 seed 1 edge 19| init_c0 = 1e-6|
|n=500 seed 1 edge 20| init_c0 = 1e-5|
|n=500 seed 1 edge 27| init_c0 = 1e-6|

Note that the set of parameters to obtain convergence is not unique, but there are the ones that are used to obtain results in for our work.

Note that bash scripts to run jobs in bulk are also provided for convenience.

