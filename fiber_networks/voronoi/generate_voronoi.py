import numpy as np
from scipy.spatial import Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt
from copy import deepcopy
import networkx as nx
import pathlib
import sys
#----------Import functions---------------------#
import sys
sys.path.append('../../utils')
import meshing 
import graph

'''
Every function returns the mesh and graph in the respective folders
'''

def generate_network(W, H, n_points, random_state, num_segments, mesh_threshold, root_dir = '.'):
    '''
    Create a mesh of random simple chain and stores the initial coordinates of links
    Args:
        W: width of window
        H: height of window
        n_points: number Voronoi seeds
        random_state: seed for voronoi inital seed
        num_segments: number of elements per beam
        mesh_threshold: threshold fiber length to use only 2 nodes
        root_dir: root directory to save mesh
    '''
    eps = sys.float_info.epsilon
    #-------------Set seed---------------------#
    rng = np.random.RandomState(random_state)

    #------Get Line centers and orientation----#
    centers_x = rng.uniform(0,W,n_points)
    centers_y = rng.uniform(0,H,n_points)

    points = np.array([centers_x,centers_y]).T

    vor = Voronoi(points)

    vertices = deepcopy(vor.vertices)
    ridge_vertices = deepcopy(vor.ridge_vertices)

    # for polygons with points at infinity
    center = points.mean(axis=0)
    ptp_bound = vor.points.ptp(axis=0)
    # got this from source code
    for count, (pointidx, idx) in enumerate(zip(vor.ridge_points, vor.ridge_vertices)):
        idx = np.asarray(idx)
        non_inf = False
        if np.any(idx < 0): # infinite points
            ii = idx[1] # finite end Voronoi vertex
            t = points[pointidx[1]] - points[pointidx[0]]  # tangent
            t = t / np.linalg.norm(t)
            n = np.array([-t[1], t[0]]) # normal
            midpoint = points[pointidx].mean(axis=0)
            far_point = vor.vertices[ii] + np.sign(np.dot(midpoint - center, n)) * n * ptp_bound.max()

            vertices = np.append(vertices, [far_point], axis = 0)
            ridge_vertices[count][0] = len(vertices)-1 # replace the end at infinity
            end_1 = deepcopy(far_point)
            end_2 = deepcopy(vertices[ii])
        else: # finite points
            non_inf = True
            end_1 = deepcopy(vertices[idx[0]])
            end_2 = deepcopy(vertices[idx[1]])
    
                
        # Bound y-direction:
        if end_1[1] < 0 and ridge_vertices[count] != [False, False]:
            end_1[0] = end_2[0] + (-end_2[1])*((end_1[0]-end_2[0])/(end_1[1]-end_2[1]+eps))
            end_1[1] = 0.0
            # if non_inf:
            vertices = np.append(vertices, [end_1], axis = 0)
            ridge_vertices[count][0] = len(vertices)-1
        elif end_1[1] > H and ridge_vertices[count] != [False, False]:
            end_1[0] = end_2[0] + (H-end_2[1])*((end_1[0]-end_2[0])/(end_1[1]-end_2[1]+eps))
            end_1[1] = H
            # if non_inf:
            vertices = np.append(vertices, [end_1], axis = 0)
            ridge_vertices[count][0] = len(vertices)-1
        if end_2[1] < 0 and ridge_vertices[count] != [False, False]:
            end_2[0] = end_1[0] + (-end_1[1])*((end_2[0]-end_1[0])/(end_2[1]-end_1[1]+eps))
            end_2[1] = 0.0 
            # if non_inf:
            vertices = np.append(vertices, [end_2], axis = 0)
            ridge_vertices[count][1] = len(vertices)-1
        elif end_2[1] > H and ridge_vertices[count] != [False, False]:
            end_2[0] = end_1[0] + (H-end_1[1])*((end_2[0]-end_1[0])/(end_2[1]-end_1[1]+eps))
            end_2[1] = H   
            # if non_inf:
            vertices = np.append(vertices, [end_2], axis = 0)
            ridge_vertices[count][1] = len(vertices)-1 
        
        # Remove dangling edges in the x direction
        if end_1[0] < 0 :
            ridge_vertices[count] = [False, False]#np.zeros_like(ridge_vertices[count], dtype=bool)
        elif end_1[0] > W :
            ridge_vertices[count] = [False, False]#np.zeros_like(ridge_vertices[count], dtype=bool)

        if end_2[0] < 0 :
            ridge_vertices[count] = [False, False]#np.zeros_like(ridge_vertices[count], dtype=bool)
        elif end_2[0] > W :
            ridge_vertices[count] = [False, False]#np.zeros_like(ridge_vertices[count], dtype=bool)

    # Points to generate mesh
    n_ridge = len(ridge_vertices)        
    end_1_mesh = []
    end_2_mesh = []

    for ii in range(n_ridge):
            idx = np.array(ridge_vertices[ii])
            end_1_mesh.append([vertices[idx[0],0],vertices[idx[0],1]])
            end_2_mesh.append([vertices[idx[1],0],vertices[idx[1],1]])

    G = graph.create_fiber_network(vertices,ridge_vertices, prune = True)
    pos = nx.get_node_attributes(G,'pos')

    edge_dist = nx.get_edge_attributes(G,'dist').values() # get edge weights
    edge_dist = np.array(list(edge_dist))
    mean_dist = np.mean(edge_dist)

    end_1_mesh = []
    end_2_mesh = []
    for ii,jj in G.edges():
        end_1_mesh.append(pos[ii])
        end_2_mesh.append(pos[jj])

    end_1_mesh = np.array(end_1_mesh).T
    end_2_mesh = np.array(end_2_mesh).T


    mesh_points = []
    for ii in range(end_1_mesh.shape[1]):
        if (end_1_mesh[:,ii] != end_2_mesh[:,ii]).all():
            mesh_points.append([end_1_mesh[:,ii], end_2_mesh[:,ii]])
    mesh_points = np.array(mesh_points)
    # mesh
    folder = f'n{n_points}'
    f_name = f'voronoi{random_state}'
    meshing.create_network_mesh(folder,f_name, mesh_points, H, edge_dist, mean_dist,mesh_threshold, num_segments, root_dir)

    # save voronoi points to construct graph
    pathlib.Path(f'./graph/nodes/n{n_points}').mkdir(parents=True, exist_ok=True)
    pathlib.Path(f'./graph/edges/n{n_points}').mkdir(parents=True, exist_ok=True)

    np.savetxt(f'./graph/nodes/n{n_points}/seed{random_state}.txt', np.array(list(pos.values())))
    np.savetxt(f'./graph/edges/n{n_points}/seed{random_state}.txt', np.array(G.edges()), fmt = '%i')