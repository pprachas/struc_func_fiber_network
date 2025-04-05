import numpy as np
import matplotlib.pyplot as plt
import gmsh
import meshio
import pathlib

def create_mesh(folder, f_name, mesh_points, char_length):
    '''
    Create simple chain mesh given file name, locations of crosslinks and characteristic length
    Note: file name should have no extension

    Args:
        folder: folder to save mesh
        f_name: mesh name
        mesh_points: crosslink points
        char_length: characteristic length
    Returns:
        None; mesh will be saved to f_name
    '''
    # create folders to store .msh, .xdmf and .h5 files if they don;t exist
    pathlib.Path(f'./mesh_msh/{folder}').mkdir(parents=True, exist_ok=True)
    pathlib.Path(f'./mesh/{folder}').mkdir(parents=True, exist_ok=True)

    gmsh.initialize()

    
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin",char_length)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax",char_length)

    for ii in range(0,len(mesh_points)): #len(mesh_points)
        for jj in range(0,len(mesh_points[ii])-1):
            p1 = gmsh.model.occ.addPoint(mesh_points[ii][jj,0],mesh_points[ii][jj,1],0)
            p2 = gmsh.model.occ.addPoint(mesh_points[ii][jj+1,0],mesh_points[ii][jj+1,1],0)
            line=gmsh.model.occ.addLine(p1,p2)
            if ii != 0 :
                gmsh.model.occ.fragment([gmsh.model.occ.getEntities(1)[0]],gmsh.model.occ.getEntities(1)[1:])

    gmsh.model.occ.synchronize()

    gmsh.model.occ.fragment([gmsh.model.occ.getEntities(1)[0]],gmsh.model.occ.getEntities(1)[1:])

    network = gmsh.model.occ.getEntities(1)

    print(len(network))
    for ii in range(len(network)):
        group = network[ii][1]
        gmsh.model.addPhysicalGroup(1, [group], ii) # make sure each fiber is one group

    msh = gmsh.model.mesh.generate(dim=1)
    gmsh.write(f'mesh_msh/{folder}/{f_name}.msh')

    mesh = meshio.read(f'mesh_msh/{folder}/{f_name}.msh')

    mesh.points = mesh.points[:, :2] #prune z = 0 for 2D mesh
        
    line_cells = np.array([None])

    for cell in mesh.cells:
        if cell.type == "line":
            if line_cells.all() == None:
                line_cells = cell.data
            else:
                line_cells = np.concatenate((line_cells,cell.data))

    for key in mesh.cell_data_dict['gmsh:physical'].keys():
        if key == 'line':
            line_data = mesh.cell_data_dict['gmsh:physical'][key]        
    fiber = meshio.Mesh(points=mesh.points, cells={"line": line_cells}, cell_data={'fibers':[line_data]})
    meshio.write(f'mesh/{folder}/{f_name}.xdmf', fiber)
    gmsh.finalize()


def create_network_mesh(folder, f_name, mesh_points, char_length, edge_dist,mean_dist, mesh_threshold, num_segments,root_dir = '.'):
    '''
    Create fiber network mesh given file name, locations of crosslinks and characteristic length
    Note: file name should have no extension. Come with more refined meshing strategy compared to
    the function create_mesh

    Args:
        folder: folder to save mesh
        f_name: mesh name
        mesh_points: crosslink points
        char_length: characteristic length
        edge_dist: list of edge distances (fiber length)
        mean_dist: mean of edge dist
        mesh_threshold: threshold to mesh fiber with only 2 nodes
        num_segments: number of segments per fiber
        root_dir: root directory to store mesh
    Returns:
        None; mesh will be saved to f_name

    '''
    # create folders to store .msh, .xdmf and .h5 files if they don;t exist
    pathlib.Path(f'{root_dir}/mesh_msh/{folder}').mkdir(parents=True, exist_ok=True)
    pathlib.Path(f'{root_dir}/mesh/{folder}').mkdir(parents=True, exist_ok=True)

    gmsh.initialize()
    
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin",char_length)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax",char_length)

    lines = []

    for ii in range(0,len(mesh_points)): #len(mesh_points)
        p1 = gmsh.model.occ.addPoint(mesh_points[ii][0,0],mesh_points[ii][0,1],0, char_length)
        p2 = gmsh.model.occ.addPoint(mesh_points[ii][1,0],mesh_points[ii][1,1],0, char_length)
        line=gmsh.model.occ.addLine(p1,p2)
        # lines.append(line)

        points = []

        if edge_dist[ii]>mesh_threshold*mean_dist:
            len_segment = np.linalg.norm(mesh_points[ii][1] - mesh_points[ii][0])/num_segments
            direction = (mesh_points[ii][1] - mesh_points[ii][0])/np.linalg.norm(mesh_points[ii][1] - mesh_points[ii][0])
            segment_steps = np.arange(0,num_segments+1)
            init_point = mesh_points[ii][0]
            for kk in range(num_segments-1):
                point = init_point + segment_steps[kk+1]*len_segment*direction
                p = gmsh.model.occ.addPoint(point[0],point[1],0)
                points.append(p)
            
            # print(line)
            combined_line = gmsh.model.occ.fragment(list(zip([0]*len(points),points)),[(1,line)])[0]

            idx = [ll for ll, (v, *_) in enumerate(combined_line) if v == 1]
            # print(combined_line)
            # print([combined_line[mm][1] for mm in idx])

            
            lines.append([combined_line[mm][1] for mm in idx])

        else:
            lines.append(line)
            pass

    if ii != 1:
        gmsh.model.occ.fragment([gmsh.model.occ.getEntities(1)[0]],gmsh.model.occ.getEntities(1)[1:])

    gmsh.model.occ.fragment([gmsh.model.occ.getEntities(1)[0]],gmsh.model.occ.getEntities(1)[1:])
    gmsh.model.occ.synchronize()

    for ii in range(len(lines)):
        gmsh.model.addPhysicalGroup(1, lines[ii], ii) # make sure each fiber is one group

    msh = gmsh.model.mesh.generate(dim=1)
    gmsh.write(f'{root_dir}/mesh_msh/{folder}/{f_name}.msh')

    mesh = meshio.read(f'{root_dir}/mesh_msh/{folder}/{f_name}.msh')

    mesh.points = mesh.points[:, :2] #prune z = 0 for 2D mesh
        
    line_cells = np.array([None])

    for cell in mesh.cells:
        if cell.type == "line":
            if line_cells.all() == None:
                line_cells = cell.data
            else:
                line_cells = np.concatenate((line_cells,cell.data))

    for key in mesh.cell_data_dict['gmsh:physical'].keys():
        if key == 'line':
            line_data = mesh.cell_data_dict['gmsh:physical'][key]        
    fiber = meshio.Mesh(points=mesh.points, cells={"line": line_cells}, cell_data={'fibers':[line_data]})
    meshio.write(f'{root_dir}/mesh/{folder}/{f_name}.xdmf', fiber)
    gmsh.finalize()


def create_ablated_mesh(folder, f_name, og_mesh_name, n, seed, edges, root_dir = '.'):
    '''
    Create ablated fiber network, the original fiber network mesh must already exist

    Args:
        folder: folder to save mesh
        f_name: mesh name
        mesh_name: name of original mesh
        n: n Voronoi seeds
        seed: random seed to generate voronoi diagram
        edges: edges to ablate (as a list)

    Returns:
        None: mesh will be saved to f_name
    '''
    #----------Create folder to save mesh--------------------#
    pathlib.Path(f'{root_dir}/mesh_msh/{folder}').mkdir(parents=True, exist_ok=True)
    pathlib.Path(f'{root_dir}/mesh/{folder}').mkdir(parents=True, exist_ok=True)

    #---------- Import mesh---------------#

    gmsh.initialize()
    gmsh.open(og_mesh_name)

    remove_edge = [(1,ii) for ii in edges]
    print(remove_edge)

    gmsh.model.removePhysicalGroups(remove_edge) # remove edge

    msh = gmsh.model.mesh.generate(dim=1)
    gmsh.write(f'{root_dir}/mesh_msh/{folder}/{f_name}.msh')

    mesh = meshio.read(f'{root_dir}/mesh_msh/{folder}/{f_name}.msh')

    mesh.points = mesh.points[:, :2] #prune z = 0 for 2D mesh
        
    line_cells = np.array([None])

    for cell in mesh.cells:
        if cell.type == "line":
            if line_cells.all() == None:
                line_cells = cell.data
            else:
                line_cells = np.concatenate((line_cells,cell.data))

    for key in mesh.cell_data_dict['gmsh:physical'].keys():
        if key == 'line':
            line_data = mesh.cell_data_dict['gmsh:physical'][key]        
    fiber = meshio.Mesh(points=mesh.points, cells={"line": line_cells}, cell_data={'fibers':[line_data]})
    meshio.write(f'{root_dir}/mesh/{folder}/{f_name}.xdmf', fiber)

    gmsh.finalize()