#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plot
import sys

data = np.genfromtxt("kapla.dat",delimiter=' ',names=['N','kapla'])

# Grafico N x kapla
fig = plot.figure()
axis = fig.add_subplot(111)
axis.set_title("N x kapla")
axis.set_xlabel('N')
axis.set_ylabel('kapla')
axis.plot(data['N'],data['kapla'],c='r',label='kapla')
leg = axis.legend()
plot.grid()
plot.show()
