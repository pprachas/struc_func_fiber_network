import numpy as np
import pathlib
import matplotlib.pyplot as plt
import networkx as nx
import sys
import time

#--------Import functions---------------------#
import sys
sys.path.append('../../utils')
sys.path.append('..')

import fea_HR as fea
import graph
from generate_voronoi import generate_network

t0 = time.time()
n_voronoi = [100,200,300,400,500] # number of seed
L = 10000 # total length of fiber
W=H=L
mesh_threshold= 0 
num_segments = 20

crit_strain_all = []
init_contour = []
num_run = 20

for n in n_voronoi:
    for seed in range(num_run):
# Name of mesh
        mesh_name= f'mesh/n{n}/voronoi{seed}.xdmf'

        # file names to write results
        f_crit_name = f'phase_diagram/crit_strain/n{int(n)}.txt'


        # run mesh
        generate_network(W,H,n,seed,num_segments,mesh_threshold)
