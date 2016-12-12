from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from matplotlib.mlab import griddata
import matplotlib.pyplot as plt
import numpy as np
import sys

data = np.genfromtxt("solution.dat")
x = data[:,0]
y = data[:,1]
z = data[:,2]
xi = np.unique(x)
yi = np.unique(y)
X, Y = np.meshgrid(xi,yi)
Z = griddata(x,y,z,xi,yi,interp='linear')

fig = plt.figure()
ax = fig.gca(projection='3d')
# Passa X, Y, Z e definicoes como tamanho do ponto, cor da legenda
surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.jet,
                    linewidth=0, antialiased=False)

plt.show()

'''
# -----------------------------------------------------------------------------------------
# Ler o arquivo que acabamos de escrever e plotar em um grafico 3d
file_reader = open('data.dat','r')

# Lendo cada linha do arquivo
for line in file_reader:
    # Quebrar a linha para cada espaco que for encontrado --> Conseguimos (x,y,z)    
    token = line.split(' ')
    # Tratar os dados ...
    
file_reader.close()
# -----------------------------------------------------------------------------------------

# -----------------------------------------------------------------------------------------
# Plotar o grafico
fig = plt.figure()
ax = fig.gca(projection='3d')
# Passa X, Y, Z e definicoes como tamanho do ponto, cor da legenda
surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

plt.show()
# -----------------------------------------------------------------------------------------
'''

# Este eh um exemplo de comentario em bloco no Python
''''
fig = plt.figure()
ax = fig.gca(projection='3d')
X = np.arange(-5, 5, 0.25)
Y = np.arange(-5, 5, 0.25)
X, Y = np.meshgrid(X, Y)
R = np.sqrt(X**2 + Y**2)
Z = np.sin(R)
surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
ax.set_zlim(-1.01, 1.01)

ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()

'''