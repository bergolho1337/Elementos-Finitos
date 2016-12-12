#!/usr/bin/python

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
import numpy as np
import sys

dt = float(sys.argv[1])
tMax = float(sys.argv[2])
N = int(tMax/dt)

for k in range(N):
  if (k % 10 == 0):
    filename = "solution"+str(k)+".dat"
    data = np.genfromtxt(filename,delimiter=' ',names=['x','analit','aprox'])

    # Grafico da solucao analitica vs aproximacao
    fig = plt.figure()
    axis = fig.gca()
    axis.set_title("Equacao do calor (N = "+str(N)+")")
    axis.set_xlabel('x')
    axis.set_ylabel('u')
    axis.set_ylim(-1.01, 1.01)
    axis.plot(data['x'],data['analit'],c='r',label='analit')
    axis.plot(data['x'],data['aprox'],'.',c='b',label='aprox')
    leg = axis.legend(loc=2)
    plt.grid()
    plt.savefig('/home/lucas/Documentos/Mestrado/Metodo_Elementos_Finitos/Equacao_Calor_Transiente/Graphics/solution'+str(k)+".png")      

#plt.show()
