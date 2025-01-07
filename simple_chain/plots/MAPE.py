import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, InsetPosition, mark_inset, zoomed_inset_axes
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import pathlib
from scipy.stats import entropy


#--------Import functions---------------------#
import sys
sys.path.append('../../utils')
import graph


#-----------------Triangular Chain------------------------------#
f_name = './MAPE_data/triangular_chain/'
pathlib.Path(f_name).mkdir(parents=True, exist_ok=True)
# plot settings
plt.style.use('jeff_style.mplstyle')


kappa_tildes = [1e-3,1e-4,1e-5,1e-6]

L_0 = 10000 # total length of fiber
lmbdas = L_0/np.array([5,10,20,40]) # wavelength 
amplitudes = [L_0/10,L_0/20,L_0/40] #amplitude
E = 1.0 # Stiffness moduli

# marker = ['s','D','o','*']

# plt.figure()
# plt.xlabel(r'$\theta_0$')
# plt.ylabel(r'MAPE')
# plt.title('Triangular Chains')
# Load data

error_stretch_all = []
error_stretch_bend_all = []

area_stretch_all = []
area_stretch_bend_all = []

for count,kappa_tilde in enumerate(kappa_tildes):
    for lmbda in lmbdas:
        for num_run,amplitude in enumerate(amplitudes):


            link_points = np.loadtxt(f'../triangular_chain/link_points/a{int(int(amplitude))}/lmbda{int(lmbda)}.txt')
            crit_strain = np.loadtxt(f'../triangular_chain/phase_diagram/crit_strain/a{int(amplitude)}.txt')[num_run][count]

            G = graph.create_chain_graph(link_points)

            edge_dist = nx.get_edge_attributes(G,'dist').values() # get edge weights
            edge_dist = np.array(list(edge_dist))
            init_contour = np.sum(edge_dist)
            L_char = (lmbda/(2*L_0))*init_contour
            r = 2*L_char*np.sqrt(kappa_tilde)
            S = np.pi*r**2
            I = np.pi*r**4/4
            f_prefix_name = f'../triangular_chain/fea_run/a{int(amplitude)}/kappa_tilde{kappa_tilde}/lmbda{int(lmbda)}'
            disp = np.loadtxt(f'{f_prefix_name}/disp.txt')
            cont = np.loadtxt(f'{f_prefix_name}/cont_total.txt')
            bend = np.loadtxt(f'{f_prefix_name}/bend_total.txt')
            stretch = np.loadtxt(f'{f_prefix_name}/stretch_total.txt')
            shear = np.loadtxt(f'{f_prefix_name}/shear_total.txt')
            force = np.loadtxt(f'{f_prefix_name}/force.txt')

            strain=disp/L_0
            straightness=cont/L_0 #straightness.append((L+disp)/cont)
            bend_energy=bend
            stretch_energy=stretch
            shear_energy=shear
            total_energy=bend+stretch+shear
            normalized_force=force/(E*S)

            init_straightness = init_contour/L_0
            init_theta = np.arccos(init_straightness**(-1))

            # stretching component
            stretching = (straightness/init_straightness)-1
            L = L_0*(1+strain)

            #bending component
            sin_theta = np.sqrt((straightness*L_0)**2-L**2)/(straightness*L_0)
            theta = np.arccos(L/(straightness*L_0))
            bending = 3*kappa_tilde*(init_theta-theta)*sin_theta

            error_stretch = np.mean(np.abs((normalized_force - stretching)/normalized_force))

            error_bend_stretch = np.mean(np.abs((normalized_force - (stretching+bending))/normalized_force))

            print(np.abs(normalized_force - stretching))
            print(strain/crit_strain)

            area_stretch = np.trapz(np.abs(normalized_force - stretching), strain/crit_strain)
            area_stretch_bend = np.trapz(np.abs(normalized_force - (stretching+bending)), strain/crit_strain)

            error_stretch_all.append([init_theta,error_stretch])
            error_stretch_bend_all.append([init_theta,error_bend_stretch])

            area_stretch_all.append([init_theta, area_stretch])
            area_stretch_bend_all.append([init_theta, area_stretch_bend])

np.savetxt(f'{f_name}error_stretch.txt',np.array(error_stretch_all))
np.savetxt(f'{f_name}error_stretch_bend.txt',np.array(error_stretch_bend_all))

np.savetxt(f'{f_name}area_stretch.txt',np.array(area_stretch_all))
np.savetxt(f'{f_name}area_stretch_bend.txt',np.array(area_stretch_bend_all))

#-----------------Sinusoidal Chain------------------------------#
f_name = './MAPE_data/sinusoidal_chain/'
pathlib.Path(f_name).mkdir(parents=True, exist_ok=True)
# plot settings

kappa_tildes = [1e-3,1e-4,1e-5,1e-6]

L_0 = 10000 # total length of fiber
lmbdas = L_0/np.array([5,10,20,40]) # wavelength 
amplitudes = [L_0/10,L_0/20,L_0/40] #amplitude
E = 1.0 # Stiffness moduli
# plt.figure()
# plt.xlabel(r'$\theta_0$')
# plt.ylabel(r'MAPE')
# plt.title('Sinusoidal Chains')

error_stretch_all = []
error_stretch_bend_all = []

area_stretch_all = []
area_stretch_bend_all = []

# Load data
for count,kappa_tilde in enumerate(kappa_tildes):
    for lmbda in lmbdas:
        for num_run,amplitude in enumerate(amplitudes):


            link_points = np.loadtxt(f'../sinusoidal_chain/link_points/a{int(int(amplitude))}/lmbda{int(lmbda)}.txt')
            crit_strain = np.loadtxt(f'../triangular_chain/phase_diagram/crit_strain/a{int(amplitude)}.txt')[num_run][count]

            G = graph.create_chain_graph(link_points)

            edge_dist = nx.get_edge_attributes(G,'dist').values() # get edge weights
            edge_dist = np.array(list(edge_dist))
            init_contour = np.sum(edge_dist)
            L_char = (lmbda/(2*L_0))*init_contour
            r = 2*L_char*np.sqrt(kappa_tilde)
            S = np.pi*r**2
            I = np.pi*r**4/4
            f_prefix_name = f'../sinusoidal_chain/fea_run/a{int(amplitude)}/kappa_tilde{kappa_tilde}/lmbda{int(lmbda)}'
            disp = np.loadtxt(f'{f_prefix_name}/disp.txt')
            cont = np.loadtxt(f'{f_prefix_name}/cont_total.txt')
            bend = np.loadtxt(f'{f_prefix_name}/bend_total.txt')
            stretch = np.loadtxt(f'{f_prefix_name}/stretch_total.txt')
            shear = np.loadtxt(f'{f_prefix_name}/shear_total.txt')
            force = np.loadtxt(f'{f_prefix_name}/force.txt')

            strain=disp/L_0
            straightness=cont/L_0 #straightness.append((L+disp)/cont)
            bend_energy=bend
            stretch_energy=stretch
            shear_energy=shear
            total_energy=bend+stretch+shear
            normalized_force=force/(E*S)

            init_straightness = init_contour/L_0
            init_theta = np.arccos(init_straightness**(-1))

            # stretching component
            stretching = (straightness/init_straightness)-1
            L = L_0*(1+strain)

            #bending component
            sin_theta = np.sqrt((straightness*L_0)**2-L**2)/(straightness*L_0)
            theta = np.arccos(L/(straightness*L_0))
            bending = 3*kappa_tilde*(init_theta-theta)*sin_theta

            error_stretch = np.mean(np.abs((normalized_force - stretching)/normalized_force))

            error_bend_stretch = np.mean(np.abs((normalized_force - (stretching+bending))/normalized_force))

            area_stretch = np.trapz(np.abs(normalized_force - stretching), strain/crit_strain)
            area_stretch_bend = np.trapz(np.abs(normalized_force - (stretching+bending)), strain/crit_strain)

            error_stretch_all.append([init_theta,error_stretch])
            error_stretch_bend_all.append([init_theta,error_bend_stretch])

            area_stretch_all.append([init_theta, area_stretch])
            area_stretch_bend_all.append([init_theta, area_stretch_bend])

print(area_stretch_all)
np.savetxt(f'{f_name}error_stretch.txt',np.array(error_stretch_all))
np.savetxt(f'{f_name}error_stretch_bend.txt',np.array(error_stretch_bend_all))

np.savetxt(f'{f_name}area_stretch.txt',np.array(area_stretch_all))
np.savetxt(f'{f_name}area_stretch_bend.txt',np.array(area_stretch_bend_all))

#---------------Random Chain-----------------------------------#
f_name = './MAPE_data/random_chain/'
pathlib.Path(f_name).mkdir(parents=True, exist_ok=True)

num_links = [10,20,30,40,50,60, 70, 80] # number of links

widths = [L_0/10,L_0/20,L_0/40,L_0/50]

num_runs = 20

error_stretch_all = []
error_stretch_bend_all = []

area_stretch_all = []
area_stretch_bend_all = []

KL_all = []

#plt.figure()
# plt.xlabel(r'$\theta_0$')
# plt.ylabel(r'MAPE')
# plt.title('Random Chains')
# fig,ax = plt.subplots(1,1)
for count,kappa_tilde in enumerate(kappa_tildes):
    for num_link in num_links:
        for w in widths:
            for num_run in range(num_runs):

                link_points = np.loadtxt(f'../random_chain/link_points/w{int(w)}/n{num_link}/link_points{num_run}.txt')
                crit_strain = np.loadtxt(f'../random_chain/phase_diagram/crit_strain/w{int(w)}/n{num_link}.txt')[num_run][count]

                G = graph.create_chain_graph(link_points)
                orientation = np.fromiter(nx.get_edge_attributes(G,'orientation').values(), dtype=float)

                edge_dist = nx.get_edge_attributes(G,'dist').values() # get edge weights
                edge_dist = np.array(list(edge_dist))
                init_contour = np.sum(edge_dist)
                L_char = np.mean(edge_dist)
                r = 2*L_char*np.sqrt(kappa_tilde)
                S = np.pi*r**2
                I = np.pi*r**4/4
                f_prefix_name = f'../random_chain/fea_run/w{int(w)}/n{num_link}/kappa_tilde{kappa_tilde}/seed{int(num_run)}'
                disp = np.loadtxt(f'{f_prefix_name}/disp.txt')
                cont = np.loadtxt(f'{f_prefix_name}/cont_total.txt')
                bend = np.loadtxt(f'{f_prefix_name}/bend_total.txt')
                stretch = np.loadtxt(f'{f_prefix_name}/stretch_total.txt')
                shear = np.loadtxt(f'{f_prefix_name}/shear_total.txt')
                force = np.loadtxt(f'{f_prefix_name}/force.txt')

                strain=disp/L_0
                straightness=cont/L_0 #straightness.append((L+disp)/cont)
                bend_energy=bend
                stretch_energy=stretch
                shear_energy=shear
                total_energy=bend+stretch+shear
                normalized_force=force/(E*S)

                init_straightness = init_contour/L_0
                init_theta = np.arccos(init_straightness**(-1))


                # stretching component
                stretching = (straightness/init_straightness)-1
                L = L_0*(1+strain)

                #bending component
                sin_theta = np.sqrt((straightness*L_0)**2-L**2)/(straightness*L_0)
                theta = np.arccos(L/(straightness*L_0))
                bending = 3*kappa_tilde*(init_theta-theta)*sin_theta

                error_stretch = np.mean(np.abs((normalized_force - stretching)/normalized_force))
                error_bend_stretch = np.mean(np.abs((normalized_force - (stretching+bending))/normalized_force))

                area_stretch = np.trapz(np.abs(normalized_force - stretching), strain/crit_strain)
                area_stretch_bend = np.trapz(np.abs(normalized_force - (stretching+bending)), strain/crit_strain)

                # Symmetry of orientation
                trans1 = (orientation-np.pi/2) % (2*np.pi) # make sure angles are [0,2pi]
                trans2 = (-orientation+np.pi/2) % (2*np.pi) # make sure angles are [0,2pi]
                symm = entropy(np.sort(trans1), np.sort(trans2))
                error_stretch_all.append([init_theta,error_stretch])
                error_stretch_bend_all.append([init_theta,error_bend_stretch])
                KL_all.append(entropy(np.sort(trans1), np.sort(trans2)))

                area_stretch_all.append([init_theta, area_stretch])
                area_stretch_bend_all.append([init_theta, area_stretch_bend])

            

np.savetxt(f'{f_name}error_stretch.txt',np.array(error_stretch_all))
np.savetxt(f'{f_name}error_stretch_bend.txt',np.array(error_stretch_bend_all))
np.savetxt(f'{f_name}KL.txt', np.array(KL_all))

np.savetxt(f'{f_name}area_stretch.txt',np.array(area_stretch_all))
np.savetxt(f'{f_name}area_stretch_bend.txt',np.array(area_stretch_bend_all))

#---------Discretized Sine------------------------#
f_name = './MAPE_data/discretized_sin/'
pathlib.Path(f_name).mkdir(parents=True, exist_ok=True)
# plot settings
plt.style.use('jeff_style.mplstyle')

kappa_tildes = [1e-3,1e-4,1e-5,1e-6]

L_0 = 10000 # total length of fiber
lmbdas = L_0/np.array([1,5,10,20,40]) # wavelength 
amplitudes = [L_0/10,L_0/20,L_0/40] #amplitude
E = 1.0 # Stiffness moduli
# plt.figure()
# plt.xlabel(r'$\theta_0$')
# plt.ylabel(r'MAPE')
# plt.title('Sinusoidal Chains')

error_stretch_all = []
error_stretch_bend_all = []

area_stretch_all = []
area_stretch_bend_all = []

# Load data
for count,kappa_tilde in enumerate(kappa_tildes):
    for lmbda in lmbdas:
        for num_run,amplitude in enumerate(amplitudes):

            link_points = np.loadtxt(f'../discretized_sin/link_points/a{int(int(amplitude))}/lmbda{int(lmbda)}.txt')
            crit_strain = np.loadtxt(f'../discretized_sin/phase_diagram/crit_strain/a{int(amplitude)}.txt')[num_run][count]
            

            G = graph.create_chain_graph(link_points)

            edge_dist = nx.get_edge_attributes(G,'dist').values() # get edge weights
            edge_dist = np.array(list(edge_dist))
            init_contour = np.sum(edge_dist)
            L_char = np.mean(edge_dist)
            r = 2*L_char*np.sqrt(kappa_tilde)
            S = np.pi*r**2
            I = np.pi*r**4/4
            f_prefix_name = f'../discretized_sin/fea_run/a{int(amplitude)}/kappa_tilde{kappa_tilde}/lmbda{int(lmbda)}'
            disp = np.loadtxt(f'{f_prefix_name}/disp.txt')
            cont = np.loadtxt(f'{f_prefix_name}/cont_total.txt')
            bend = np.loadtxt(f'{f_prefix_name}/bend_total.txt')
            stretch = np.loadtxt(f'{f_prefix_name}/stretch_total.txt')
            shear = np.loadtxt(f'{f_prefix_name}/shear_total.txt')
            force = np.loadtxt(f'{f_prefix_name}/force.txt')

            strain=disp/L_0
            straightness=cont/L_0 #straightness.append((L+disp)/cont)
            bend_energy=bend
            stretch_energy=stretch
            shear_energy=shear
            total_energy=bend+stretch+shear
            normalized_force=force/(E*S)

            init_straightness = init_contour/L_0
            init_theta = np.arccos(init_straightness**(-1))

            # stretching component
            stretching = (straightness/init_straightness)-1
            L = L_0*(1+strain)

            #bending component
            sin_theta = np.sqrt((straightness*L_0)**2-L**2)/(straightness*L_0)
            theta = np.arccos(L/(straightness*L_0))
            bending = 3*kappa_tilde*(init_theta-theta)*sin_theta

            error_stretch = np.mean(np.abs((normalized_force - stretching)/normalized_force))

            error_bend_stretch = np.mean(np.abs((normalized_force - (stretching+bending))/normalized_force))

            area_stretch = np.trapz(np.abs(normalized_force - stretching), strain/crit_strain)
            area_stretch_bend = np.trapz(np.abs(normalized_force - (stretching+bending)), strain/crit_strain)

            error_stretch_all.append([init_theta,error_stretch])
            error_stretch_bend_all.append([init_theta,error_bend_stretch])

            area_stretch_all.append([init_theta, area_stretch])
            area_stretch_bend_all.append([init_theta, area_stretch_bend])


np.savetxt(f'{f_name}error_stretch.txt',np.array(error_stretch_all))
np.savetxt(f'{f_name}error_stretch_bend.txt',np.array(error_stretch_bend_all))

np.savetxt(f'{f_name}area_stretch.txt',np.array(area_stretch_all))
np.savetxt(f'{f_name}area_stretch_bend.txt',np.array(area_stretch_bend_all))