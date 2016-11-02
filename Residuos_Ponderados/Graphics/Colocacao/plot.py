import numpy as np
import matplotlib.pyplot as plot
import sys

plot.plot(*np.loadtxt("f.dat",unpack=True), linewidth=2.0, label='f')
plot.plot(*np.loadtxt("f_aprox.dat",unpack=True), linewidth=2.0, label='f')
plot.legend(['u(x)','N=2'])
plot.title("MWR usando Colocacao")
plot.grid()
plot.show()

plot.plot(*np.loadtxt("residue.dat",unpack=True), linewidth=2.0, label='f')
plot.legend(['N=2'])
plot.title("MWR usando Colocacao - Residuo")
plot.grid()
plot.show()
