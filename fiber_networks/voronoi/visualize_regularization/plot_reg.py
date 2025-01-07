import numpy as np
import matplotlib.pyplot as plt

force = np.loadtxt('force.txt')
stretch_total = np.loadtxt('stretch_total.txt')
damp_total = np.loadtxt('damp_total.txt')

plt.figure()
plt.plot(damp_total)
plt.title('Damping Energy')
plt.xlabel('Iteration number')
plt.ylabel('Energy')
plt.savefig('damp.png')


plt.figure()
plt.plot(force)
plt.title('Force')
plt.xlabel('Iteration number')
plt.ylabel('Force')
plt.savefig('force.png')

plt.figure()
plt.plot(stretch_total)
plt.title('stretching Energy')
plt.xlabel('Iteration number')
plt.ylabel('Energy')
plt.savefig('stretch.png')
