import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from scipy.spatial.distance import pdist, squareform
from itertools import islice

def create_chain_graph(link_points):
    '''
    Function to create graph of single fiber chains

    Args:
        link_points: coordinated of link points
    Returns: 
        G: Networkx graph of single fiber chains
    '''
    G = nx.Graph()
    
    # create vertices
    vertices_dict = []
    for ii in link_points:
        vertices_dict.append({'pos':ii})
    
    vertices_dict = list(enumerate(vertices_dict))

    # create edges
    edges_idx = np.arange(0,len(link_points))
    edges = []
    for ii in range(len(edges_idx)-1):
        edges.append([edges_idx[ii],edges_idx[ii+1]])
    
    G.add_nodes_from(vertices_dict)
    G.add_edges_from(edges)

    # add edge weights (Euclidean distance)
    pos_array = []
    for node,nodal_pos in G.nodes(data = 'pos'):  
        pos_array.append(nodal_pos)

    distances = squareform(pdist(np.array(pos_array))) # matrix form of pairwise distance
    
    for ii,jj in G.edges():
        G[ii][jj]['dist'] = distances[ii][jj]
        G[ii][jj]['orientation'] = np.arccos((G.nodes()[jj]['pos'] - G.nodes()[ii]['pos'])[0]/G[ii][jj]['dist'])

    return G

def create_fiber_network(nodes, edges, prune = False):
    '''
    Function to create graph of fiber network

    Args:
        nodes: nodes or voronoi diagram
        edges: edges of voronoi diagram
    Returns: 
        G: Networkx graph of single fiber chains
    '''
    nodes_dict = []
    for ii in nodes:
        nodes_dict.append({'pos' : ii})

    nodes_dict = list(enumerate(nodes_dict))
   
    # Initialize and create graph
    G = nx.Graph()
    G.add_nodes_from(nodes_dict)
    G.add_edges_from(edges)

    # remove self loops
    G.remove_edges_from(nx.selfloop_edges(G))

    # remove isolated nodes (no edges)
    if prune:
        G = G.subgraph(max(nx.connected_components(G), key=len))
        G = nx.convert_node_labels_to_integers(G) # make all node numbering consecutive

    pos_array = []
    for node,nodal_pos in G.nodes(data = 'pos'):
        pos_array.append(nodal_pos)

    # compute distance edge distance and add as edge weight
    distances = squareform(pdist(np.array(pos_array))) 

    for count,(ii,jj) in enumerate(G.edges()):
        G[ii][jj]['dist'] = distances[ii][jj]
        G[ii][jj]['index'] = count

    
    return G

def add_boundary_nodes(G, bottom_boundary, top_boundary, width):
    '''
    Function to add nodes for shortest paths analysis
    Args:
        G: Networkx graph of voronoi network
        bottom_boundary: y-coordinate of botton boundary
        top_boundary: y-coordinate of top boundary
        width: width of fiber network (window size)

    Returns: 
        G: Networkx graph of single fiber chains
    '''
    
    # Note - only for y axis
    bot = []
    top = []
    n = max(list(G.nodes))
    for node,val in G.nodes(data = 'pos'):
        if np.abs(val[1]- bottom_boundary) < 1e-6:
            bot.append(node)
        elif np.abs(val[1]-top_boundary) < 1e-6:
            top.append(node)

    G.add_node(n+1, pos = np.array([width/2, bottom_boundary-0.1*(top_boundary - bottom_boundary)])) # coordinates are for visualization only
    G.add_node(n+2, pos = np.array([width/2, top_boundary+0.1*(top_boundary - bottom_boundary)])) # coordinates are for visualization only
    
    for ii in bot:
        G.add_edge(n+1,ii,dist = 0) # edge-weight is 0

    for ii in top:
        G.add_edge(n+2,ii,dist = 0)  # edge-weight is 0
    
    return G


def shortest_path(G, source, target, weight):
    '''
    Function for shortest apths anaylsis

        G: Networkx graph of voronoi network
        source: source node (one of the boundary nodes)
        target: target node (one of the boundary nodes)
        weight: weight for shortest paths
    Returns: 
        paths: list of nodes for all shortest paths
    '''
    
    paths = []
    while True:
        try:
            path = nx.shortest_path(G, source, target, weight)
            path.pop(0)
            path.pop(-1)
            paths.append(path)
            G.remove_nodes_from(path)
        except:
            return paths



    
    
    