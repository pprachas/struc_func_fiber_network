# plots Directory
All the code here are used to reproduce the results in Results section for single fiber chains.

# List of code files
|File Name|Description|
----------|------------
|MAPE.py|Computes and saves area between curves of FEA results and analytical reduced order model for all chains|
|compare_paths.py|Compares the energy partitions in the whole network, shortest paths, and support network|
|energy.py|Plots percent stretch, bend, shear energy ratio and percent stretch energy ratio with change in $\tilde{\kappa}$|
|fiber_recruitment.py|Plots tangenet stiffness and visualizes fiber recruitment effect|
|network_energy.py|Visualizes energy distribution on the fiber network at $1.5\varepsilon_{crit}$|
|plot_MAPE.py|Plots area between curve results|
|plot_graph.py|Plots the constructed graph of random fiber network and example of fiber network with abalted shortest path|
|plot_mesh_refinement.py|plot results of mesh refinements|
|plot_sensitivity.py|Plot FEA results from sensitivity analysis on regularization parameter|
|shortest_path.py|Computes shortest paths for all random fiber networks and saves the paths as a list of nodes in csv file|
|visualize_sensitivity.py|visualizes fiber networks with different magnitudes of regularization|
