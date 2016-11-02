import numpy as np
import matplotlib.pyplot as plot
import sys

plot.plot(*np.loadtxt("f1.dat",unpack=True), linewidth=2.0, label='f')
plot.plot(*np.loadtxt("f_aprox1.dat",unpack=True), linewidth=2.0, label='f')
plot.plot(*np.loadtxt("f_aprox2.dat",unpack=True), linewidth=2.0, label='f')
plot.plot(*np.loadtxt("f_aprox3.dat",unpack=True), linewidth=2.0, label='f')
plot.legend(['u(x)','N=1','N=2','N=3'])
plot.title("MWR using Galerkin")
plot.grid()
plot.show()

plot.plot(*np.loadtxt("residue1.dat",unpack=True), linewidth=2.0, label='f')
plot.plot(*np.loadtxt("residue2.dat",unpack=True), linewidth=2.0, label='f')
plot.plot(*np.loadtxt("residue3.dat",unpack=True), linewidth=2.0, label='f')
plot.legend(['N=1','N=2','N=3'])
plot.title("MWR using Galerkin - Residue")
plot.grid()
plot.show()
