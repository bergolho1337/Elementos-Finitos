#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plot
import sys

N = sys.argv[1]
title_name = "Funcao Chapeu"

data = np.genfromtxt("solution.dat",delimiter=' ',names=['x','analit','aprox','residue'])

# Grafico da solucao analitica vs aproximacao
fig = plot.figure()
axis = fig.add_subplot(111)
axis.set_title(str(title_name)+" (N = "+str(N)+")")
axis.set_xlabel('x')
axis.set_ylabel('y')
axis.plot(data['x'],data['analit'],c='r',label='analit')
axis.plot(data['x'],data['aprox'],c='b',label='aprox')
leg = axis.legend(loc=2)
plot.grid()
#plot.savefig('/home/lucas/Documentos/Mestrado/Metodo_Elementos_Finitos/Projecao_L2/solucao.png')
plot.show()

# Grafico do residuo
fig = plot.figure()
axis = fig.add_subplot(111)
axis.set_title(str(title_name)+" (N = "+str(N)+")")
axis.set_xlabel('x')
axis.set_ylabel('R')
axis.plot(data['x'],data['residue'],c='g',label='residue')
leg = axis.legend(loc=2)
plot.grid()
#plot.savefig('/home/lucas/Documentos/Mestrado/Metodo_Elementos_Finitos/Projecao_L2/residuo.png')
plot.show()
