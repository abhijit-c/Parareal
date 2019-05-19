import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

T = 100
dt = .0001
tt = np.arange(0, T+dt, dt)
data = np.loadtxt("Data/fine_lorenz.txt")

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(data[:,0], data[:,1], data[:,2])
ax.set_title('Fine Lorenz system')
plt.show()

T = 100
dt = .1
tt = np.arange(0, T+dt, dt)
data = np.loadtxt("Data/parareal_lorenz.txt")

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(data[:,0], data[:,1], data[:,2])
ax.set_title('Parareal Lorenz system')
plt.show()
