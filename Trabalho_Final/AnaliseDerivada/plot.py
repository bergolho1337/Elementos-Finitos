#!/usr/bin/python

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
import numpy as np
import sys

N = int(sys.argv[1])
tipoElem = int(sys.argv[2])
if (tipoElem == 1):
  titleName = "Linear"
else:
  titleName = "Hermite"

filename = "solution.dat"
data = np.genfromtxt(filename,delimiter=' ',names=['x','analit','aprox'])

# Grafico da solucao analitica vs aproximacao
fig = plt.figure()
axis = fig.gca()
axis.set_title(titleName+" N = "+str(N))
axis.set_xlabel('x')
axis.set_ylabel('u')
axis.plot(data['x'],data['analit'],c='r',label='analit')
axis.plot(data['x'],data['aprox'],'.',c='b',label='aprox')
leg = axis.legend(loc=2)
plt.grid()
plt.show()

filename = "derivative.dat"
data = np.genfromtxt(filename,delimiter=' ',names=['x','analit','aprox'])

# Grafico da solucao analitica vs aproximacao
fig = plt.figure()
axis = fig.gca()
axis.set_title(titleName+" N = "+str(N))
axis.set_xlabel('x')
axis.set_ylabel('du')
axis.plot(data['x'],data['analit'],c='r',label='analit')
axis.plot(data['x'],data['aprox'],'.',c='b',label='aprox')
leg = axis.legend(loc=2)
plt.grid()
plt.show()