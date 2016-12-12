#!/usr/bin/python

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
import numpy as np

# Grafico da solucao analitica
fig = plt.figure()
ax = fig.gca(projection='3d')

data = np.genfromtxt("analit.dat")
x = data[:,0]
y = data[:,1]
z = data[:,2]

xi = np.linspace(min(x),max(x))
yi = np.linspace(min(y),max(y))

X, Y = np.meshgrid(xi,yi)
Z = griddata(x,y,z,xi,yi)

surf = ax.plot_surface(X,Y,Z,rstride=5,cstride=5,linewidth=1,antialiased=True)

ax.set_zlim3d(np.min(Z),np.max(Z))

plt.show()

# Grafico da solucao aproximada
fig = plt.figure()
ax = fig.gca(projection='3d')

data = np.genfromtxt("solution.dat")
x = data[:,0]
y = data[:,1]
z = data[:,2]

xi = np.linspace(min(x),max(x))
yi = np.linspace(min(y),max(y))

X, Y = np.meshgrid(xi,yi)
Z = griddata(x,y,z,xi,yi)

surf = ax.plot_surface(X,Y,Z,rstride=5,cstride=5,linewidth=1,antialiased=True,color='r')

ax.set_zlim3d(np.min(Z),np.max(Z))

plt.show()
