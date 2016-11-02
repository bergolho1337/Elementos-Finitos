import numpy as np
import matplotlib.pyplot as plot
import sys

plot.plotfile('f_aprox.dat', delimiter=' ', cols=(0, 1),names=('x', 'f'), marker='o', label='f_aprox')
plot.plot(*np.loadtxt("f.dat",unpack=True), linewidth=2.0, label='f')
plot.grid()
plot.show()

plot.plotfile('residue.dat', delimiter=' ', cols=(0, 1),names=('x', 'R'), marker='o', label='residuo')
plot.grid()
plot.show()
