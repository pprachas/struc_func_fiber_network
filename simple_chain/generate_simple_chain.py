import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import pathlib
import sys
from scipy import signal
#----------Import functions---------------------#
import sys
sys.path.append('../utils')
import meshing 
import graph


'''
Every function returns the mesh and graph in the respective folders
'''
def create_random_simple_chain(n, L, w, num_segments, seed):
    '''
    Create a mesh of random simple chain and stores the initial coordinates of links
    Args: 
        n: number of links
        L: total end-to-end length of fiber
        char_length: mesh characterictic length
        seed: seed to generate chain
    '''
    #------Randomly sample fiber locations-----#
    #Set Seed
    rng = np.random.RandomState(seed)

    #Get fiber crosslinking points
    x = rng.uniform(-w,w,n-2)
    y = rng.uniform(0,L,n-2)
    y = np.sort(y)

    mesh_points = [[[0.0,0.0], [x[0],y[0]]]]

    link_points = [[0.0,0.0]]
    for ii in range(len(x)-1):
        mesh_points.append([[x[ii],y[ii]],[x[ii+1],y[ii+1]]])
        link_points.append([x[ii],y[ii]])

    mesh_points.append([[x[-1],y[-1]],[0.0,L]])
    link_points.append([x[-1],y[-1]])
    link_points.append([0.0,L])

    mesh_points = np.array(mesh_points)
    link_points = np.array(link_points)

    pathlib.Path(f'./link_points/w{int(w)}/n{n}').mkdir(parents=True, exist_ok=True)
    np.savetxt(f'link_points/w{int(w)}/n{n}/link_points{seed}.txt', link_points)

    folder = f'w{int(w)}/n{len(link_points)}'
    f_name = f'random_chain{seed}' 

    #create graph to get contour length
    G = graph.create_chain_graph(link_points)

    edge_dist = nx.get_edge_attributes(G,'dist').values() # get edge weights
    edge_dist = np.array(list(edge_dist))
    min_dist = np.min(edge_dist)
    mean_dist = np.mean(edge_dist)

    # create mesh; mesh characteristic length depends on countour length
    meshing.create_network_mesh(folder, f_name, mesh_points, L, edge_dist, mean_dist, 0.0, num_segments)


def create_sinusoidal_simple_chain(amplitude, wavelength, num_points, L):
    '''
    Create a mesh of sinusoidal simple chain and stores the initial coordinates of links
    Args: 
        n: number of links
        L: total end-to-end length of fiber
        char_length: mesh characterictic length
        seed: seed to generate chain
    '''
    #------Create chains----#
    y = np.linspace(0,L,num_points)
    x = amplitude*np.sin(y*(2*np.pi/wavelength))

    mesh_points = []

    link_points = []
    for ii in range(len(x)-1):
        mesh_points.append([[x[ii],y[ii]],[x[ii+1],y[ii+1]]])
        link_points.append([x[ii],y[ii]])

    mesh_points.append([[x[-2],y[-2]],[x[-1],y[-1]]])
    link_points.append([x[-1],y[-1]])

    mesh_points = np.array(mesh_points)

    link_points = np.array(link_points)

    pathlib.Path(f'./link_points/a{int(amplitude)}').mkdir(parents=True, exist_ok=True)
    np.savetxt(f'./link_points/a{int(amplitude)}/lmbda{int(wavelength)}.txt', link_points)

    folder = f'a{int(amplitude)}'
    f_name = f'simple_chain_lmbda{int(wavelength)}' 

    # create mesh; mesh characteristic length depends on countour length
    meshing.create_mesh(folder,f_name, mesh_points, L)

def create_triangular_simple_chain(amplitude, wavelength, L, num_segments):
    '''
    Create a mesh of triangular wave simple chain and stores the initial coordinates of links
    Args: 
        n: number of links
        L: total end-to-end length of fiber
        char_length: mesh characterictic length
        seed: seed to generate chain
    '''
    #------Randomly sample fiber locations----#
    num_points = int(2*L/wavelength) + 1
    y = np.linspace(0,L,num_points)
    x = amplitude*signal.sawtooth(y*(2*np.pi/wavelength), width = 0.5) + L/10

    mesh_points = []

    link_points = []
    for ii in range(len(x)-1):
        mesh_points.append([[x[ii],y[ii]],[x[ii+1],y[ii+1]]])
        link_points.append([x[ii],y[ii]])

    link_points.append([x[-1],y[-1]])

    mesh_points = np.array(mesh_points)
    link_points = np.array(link_points)

    pathlib.Path(f'./link_points/a{int(amplitude)}').mkdir(parents=True, exist_ok=True)
    np.savetxt(f'./link_points/a{int(amplitude)}/lmbda{int(wavelength)}.txt', link_points)

    folder = f'a{int(amplitude)}'
    f_name = f'simple_chain_lmbda{int(wavelength)}' 

    #create graph to get contour length
    G = graph.create_chain_graph(link_points)

    edge_dist = nx.get_edge_attributes(G,'dist').values() # get edge weights
    edge_dist = np.array(list(edge_dist))

    min_dist = np.min(edge_dist)
    mean_dist = np.mean(edge_dist)

    # create mesh; mesh characteristic length depends on countour length
    meshing.create_network_mesh(folder, f_name, mesh_points, L*1e6, edge_dist, mean_dist, 0.0, num_segments)

def create_discretized_sin_simple_chain(amplitude, wavelength, num_points, L, num_segments):
    '''
    Create a mesh of sinusoidal simple chain and stores the initial coordinates of links
    Args: 
        n: number of links
        L: total end-to-end length of fiber
        char_length: mesh characteristic length
        seed: seed to generate chain
    '''
    #------Create chains----#
    y = np.linspace(0,L,num_points+1)
    x = amplitude*np.sin(y*(2*np.pi/wavelength))

    mesh_points = []

    link_points = []
    for ii in range(len(x)-1):
        mesh_points.append([[x[ii],y[ii]],[x[ii+1],y[ii+1]]])
        link_points.append([x[ii],y[ii]])

    link_points.append([x[-1],y[-1]])

    mesh_points = np.array(mesh_points)

    link_points = np.array(link_points)

    #create graph to get contour length
    G = graph.create_chain_graph(link_points)

    edge_dist = nx.get_edge_attributes(G,'dist').values() # get edge weights
    edge_dist = np.array(list(edge_dist))
    min_dist = np.min(edge_dist)

    pathlib.Path(f'./link_points/a{int(amplitude)}').mkdir(parents=True, exist_ok=True)
    np.savetxt(f'./link_points/a{int(amplitude)}/lmbda{int(wavelength)}.txt', link_points)

    folder = f'a{int(amplitude)}'
    f_name = f'simple_chain_lmbda{int(wavelength)}' 

    #create graph to get contour length
    G = graph.create_chain_graph(link_points)

    edge_dist = nx.get_edge_attributes(G,'dist').values() # get edge weights
    edge_dist = np.array(list(edge_dist))

    min_dist = np.min(edge_dist)
    mean_dist = np.mean(edge_dist)

    meshing.create_network_mesh(folder, f_name, mesh_points, L*1e6, edge_dist, mean_dist, 0.0, num_segments)