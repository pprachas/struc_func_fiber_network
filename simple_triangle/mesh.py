import numpy as np
import matplotlib.pyplot as plt
#------------Import functions-------#
import sys
sys.path.append('../utils')
import meshing

#---------Mesh Geometry---------#
# Outer points:
mesh_points = [[[0.0,0.0], [-5000,5000]], [[-5000,5000], [0,10000]],[[2500.0,0.0], [7500,5000]],[[7500,5000], [2500,10000]], [[-5000,5000], [1250,5050]], [[1250,5050],[7500,5000]]]

mesh_points = np.array(mesh_points)

plt.figure()
for ii in range(len(mesh_points)):
    plt.plot([mesh_points[ii][0,0], mesh_points[ii][1,0]],[mesh_points[ii][0,1], mesh_points[ii][1,1]], c = 'k' ,lw = 2.0)


#----------Create Mesh----------#
folder = '.'
f_name = 'simple_triangle'
char_length = 10000
edge_dist = np.ones(len(mesh_points)) # for dummy
mean_dist = 0
mesh_threshold = 0.0
num_segments = 20


meshing.create_network_mesh(folder, f_name, mesh_points, char_length, edge_dist, mean_dist, mesh_threshold, num_segments)

plt.savefig('simple_triangle.pdf')
plt.show()