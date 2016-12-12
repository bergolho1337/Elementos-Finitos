#!/usr/bin/env python
#
# Programa exemplo MEF1D 
# Problema:
#  -u_xx = f(x), em [0,1]
# com
#   f(x) = 1.0
#   u(0) = u(1) = 0
#
# Solucao exata
#   u(x) = x(1-x)/2
#
# Bernardo M. Rocha
# 16/10/2015

import sys
import numpy as np
from pylab import *

def MatrizRigidez(x):
    n = len(x)
    h = x[1] - x[0]
    A = np.zeros((n,n))
    for i in range(1,n):
        A[i-1,i-1] +=  1.0/h
        A[i-1,i]   += -1.0/h
        A[i,i-1]   += -1.0/h
        A[i,i]     +=  1.0/h
    return A

def VetorCarga(x,f):
    n = len(x)
    h = x[1] - x[0]
    b = np.zeros((n))
    for i in range(1,n):
        b[i-1] += f(x[i-1])*h/2.
        b[i]   += f(x[i])*h/2.
    return b

def CondicaoContorno(A, b, k0,kl,g0,gl):
    # impoe condicoes de contorno pelo
    # metodo da "penalizacao"
    n = len(b)
    A[0,0] += k0
    A[n-1,n-1] += kl
    b[0] += k0*g0
    b[n-1] += kl*gl

def MEF(n):

    # discretizacao de I=[0,1]
    x = np.linspace(0,1,n)
    h = x[1]-x[0]

    # funcao do lado direito f(x)
    f = lambda x: 1.0

    # monta sistema
    A = MatrizRigidez(x)
    b = VetorCarga(x,f)

    # aplica condicoes de contorno
    kappa = 10e+6
    g0, gl = 0.0, 0.0
    CondicaoContorno(A,b,kappa,kappa,g0,gl)

    print('Sistema com condicoes de contorno')
    print(A)
    print(b)

    # calcula solucao
    uh = np.linalg.solve(A,b)

    print('Solucao u')
    print(uh)

    # visualiza solucao
    subplot(211)
    x = linspace(0,1,500)
    xh = linspace(0,1,n)
    ue = x*(1-x)/2.0
    plot(xh, uh, 'r-', label='aproximacao')
    plot(x, ue, 'b-', label="exata")
    ylabel('u(x)')
    grid(True)
    legend(loc='best')
    title('Solucao u(x)')

    # visualiza derivada
    subplot(212)
    due = ((1.0 - 2.0*x)/2.0)
    duh = zeros(n-1)
    plot(x, due, 'b-', label="exata")
    for e in range(n-1):
        n0, n1 = e, e+1
        u0, u1 = uh[n0], uh[n1]
        val = (u1-u0)/h
        duh[e] = val
        # para plotar por elemento
        xx = linspace(xh[n0], xh[n1], 100)
        du = zeros(100)
        du.fill(val)
        plot(xx, du, 'r-')
    xlabel('x')
    ylabel('du(x)')
    grid(True)
    title('Derivada da solucao du(x)')
    show()

    print("Derivada du(x)")
    print(duh)

if __name__ == "__main__":

    n = int(input("Digite numero de nos: "))
    MEF(n)
