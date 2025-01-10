# utils Directory

This directory contains most of the functions used to obtain results for this work. This includes meshing, FEA analysis, and meshing.

## fea_HR.py
The main fiule for our FEA analysis, this includes regualirzation (also referred to as ``damping''). Most of the functions are used to construct the finite element problem (e.g. construct function space, mesh reading). The main functions are:

|Function|Description|
---------|------------
|run_critical_strain| Finds the critical strain transition point|
|run_fea| runs the FEA problem used to perform analyis in the work|

More details inputs and outputs of these functions can be found in the code.

## graph.py
Code to generate graphs of single fiber chains and random fiber networks. 

## meshing.py
Code to mesh single fiber chains and random fiber networks with prescribed number of elements per fiber.
