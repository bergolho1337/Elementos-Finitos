# **************************************************************************************
# Calcula as matrizes de massa e de rigidez do elemento de Hermite             __
#                                                                             /  \   
# Autor: Lucas Berg                                                           |      
#                                                                             |    _|
# Last Update: 12/12/16                                                       |   |_| x
#                                                                             |
#                                                                           \_/
# **************************************************************************************                                                                               
from sympy import *
import numpy as np

# Exemplo do Sympy
#init_printing(use_unicode=False,wrap_line=False,no_global=True)
#x = Symbol('x')
#result = integrate(x**2,x)
#pprint(result)

init_printing(use_unicode=False,wrap_line=False,no_global=True)
x = Symbol('x')
h = Symbol('h')

# Funcoes de Hermite
phi_1 = 2*(x/h)**3 - 3*(x/h)**2 + 1
phi_2 = 3*(x/h)**2 - 2*(x/h)**3
phi_3 = (x/h)**3 - 2*(x/h)**2 + (x/h)
phi_4 = (x/h)**3 - (x/h)**2
# Derivadas das funcoes de Hermite
dphi_1 = diff(phi_1,x)
dphi_2 = diff(phi_2,x)
dphi_3 = diff(phi_3,x)
dphi_4 = diff(phi_4,x)

phi = [phi_1,phi_2,phi_3,phi_4]
dphi = [dphi_1,dphi_2,dphi_3,dphi_4]

print "=========================== Local Mass Matrix ==================================="
for i in range(4):
    for j in range(4):
        print "-------------------------------------------------------------------------"
        print "phi_"+str(i)+str(j)+" = "
        pprint(integrate(phi[i]*phi[j],(x,0,h)))
        print "-------------------------------------------------------------------------"
print "================================================================================="

print

print "=========================== Local Stiff Matrix ==================================="
for i in range(4):
    for j in range(4):
        print "-------------------------------------------------------------------------"
        print "dphi_"+str(i)+str(j)+" = "
        pprint(integrate(dphi[i]*dphi[j],(x,0,h)))
        print "-------------------------------------------------------------------------"
print "================================================================================="