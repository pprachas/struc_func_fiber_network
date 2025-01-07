import matplotlib.pyplot as plt
import numpy as np
from scipy import signal

# Sinusoidal Chain
amplitude=500 
L = 10000
wavelength = [1000,2000,2500] 
num_points = 1000
n = [30]

for ii in range(1):
    y = np.linspace(0,L,num_points)
    x = amplitude*np.sin(y*(2*np.pi/wavelength[ii]))
    plt.figure()
    plt.plot(x,y, lw = 1, c = 'k')
    plt.axis('equal')
    plt.savefig(f'sin_chain{wavelength[ii]}.pdf')

# Triangular Chain
    num_points = int(2*L/wavelength[ii]) + 1
    y = np.linspace(0,L,num_points)
    x = amplitude*signal.sawtooth(y*(2*np.pi/wavelength[ii]), width = 0.5) + L/10

    plt.figure()
    plt.plot(x,y, lw = 1, c = 'k')
    plt.axis('equal')
    plt.savefig(f'tri_chain{wavelength[ii]}.pdf')

# Random Chain
    rng = np.random.RandomState(4)
    w = L/10

    #Get fiber crosslinking points
    x = rng.uniform(-w,w,n[ii]-2)
    y = rng.uniform(0,L,n[ii]-2)
    y = np.sort(y)

    plt.figure()
    plt.plot(x,y, lw = 1, c = 'k')
    plt.axis('equal')
    plt.savefig(f'random_chain{ii}.pdf')
plt.show()