#!/usr/bin/env python
#
# Programa exemplo para solucao do problema 2D:
#    -Delta(u) = f(x) em \Omega
# com
#    u = 0 sobre \partial \Omega

# usando o MEF com triangulos lineares continuos por partes
#
# Bernardo M. Rocha
# 22/04/2015

import sys, os.path
import numpy as np
from scipy import sparse
from math import sqrt

def MatrizMassa(p,t):
    npts = np.size(p,1)
    ntri = np.size(t,1)
    M = np.zeros((npts,npts))
    for k in range(ntri):
        loc2glb = t[:,k]
        x = p[0, loc2glb]
        y = p[1, loc2glb]
        area = PolyArea(x,y)
        Mk = np.array([[2,1,1],
                       [1,2,1],
                       [1,1,2]])*(area/12.)
        for il in range(3):
            i = loc2glb[il]
            for jl in range(3):
                j = loc2glb[jl]
                M[i,j] = M[i,j] + Mk[il,jl]
    return M

def MatrizRigidez(p,t,a):
    npts = np.size(p,1)
    ntri = np.size(t,1)
    A = np.zeros((npts,npts))
    for k in range(ntri):
        loc2glb = t[:,k]
        x = p[0, loc2glb]
        y = p[1, loc2glb]
        area, b, c = HatGradients(x,y)
        xc = x.mean()
        yc = y.mean()
        abar = a(xc,yc)
        Ak = abar * ( np.outer(b,b.transpose()) +
                      np.outer(c,c.transpose()) ) * area

        for il in range(3):
            i = loc2glb[il]
            for jl in range(3):
                j = loc2glb[jl]
                A[i,j] = A[i,j] + Ak[il,jl]
    return A

def VetorCarga(p,t,f):
    npts = np.size(p,1)
    ntri = np.size(t,1)
    b = np.zeros(npts)
    for k in range(ntri):
        loc2glb = t[:,k]
        x = p[0, loc2glb]
        y = p[1, loc2glb]
        area = PolyArea(x,y)
        bk = np.array([[f(x[0],y[0])],
                       [f(x[1],y[1])],
                       [f(x[2],y[2])]])*(area/3.0)
        for il in range(3):
            i = loc2glb[il]
            b[i] = b[i] + bk[il]
    return b

def RobinMatrizMassa(p,e,kappa):
    npts = np.size(p,1)
    nedg = np.size(e,0)
    R = np.zeros((npts,npts))
    for ek in range(nedg):
        loc2glb = e[ek,0:2]
        bmarker = e[ek,2]
        x = p[0, loc2glb]
        y = p[1, loc2glb]
        edgelen = sqrt( (x[1]-x[0])**2 + (y[1]-y[0])**2 )
        xc = x.mean()
        yc = y.mean()
        k = kappa(xc,yc,bmarker)
        Re = (k/6.0)*( np.array([[2.0,1.0], [1.0,2.0]]) )*edgelen
        for il in range(2):
            i = loc2glb[il]
            for jl in range(2):
                j = loc2glb[jl]
                R[i,j] = R[i,j] + Re[il,jl]
        
    return R

def RobinVetorCarga(p,e,kappa,gD,gN):
    npts = np.size(p,1)
    nedg = np.size(e,0)
    r = np.zeros(npts)
    for ek in range(nedg):
        loc2glb = e[ek,0:2]
        bmarker = e[ek,2]
        x = p[0, loc2glb]
        y = p[1, loc2glb]
        elen = sqrt((x[0]-x[1])**2 + (y[0]-y[1])**2)
        xc = x.mean()
        yc = y.mean()
        tmp = kappa(xc,yc,bmarker)*gD(xc,yc,bmarker) + gN(xc,yc,bmarker)
        re = tmp * np.array([[1],[1]])*(elen/2)
        for il in range(2):
            i = loc2glb[il]
            r[i] = r[i] + re[il]
    return r

def HatGradients(x,y):
    area = PolyArea(x,y)
    b = np.array([y[1]-y[2],
                  y[2]-y[0],
                  y[0]-y[1]])*(1./(2.*area))
    c = np.array([x[2]-x[1],
                  x[0]-x[2],
                  x[1]-x[0]])*(1./(2.*area))
    return area,b,c
        
def PolyArea(x, y):
    """
    Calcula a area de um poligono com vertices em (x,y).
    Ref: http://mathworld.wolfram.com/PolygonArea.html
    """
    area = 0.0
    npts = max(len(x), len(y))
    for ii in range(npts):
        area += x[ii]*y[(ii+1) % npts] - x[(ii+1) % npts]*y[ii]
    return  np.abs(area*0.5)

def LeMalha(basename):
    """
    Le os arquivos .node e .ele. Se disponivel le o .edge
    """
    pfile = basename + ".node"
    tfile = basename + ".ele"
    efile = basename + ".edge"
    
    pts = np.loadtxt(pfile,skiprows=1,comments='#')
    pts = pts[:,1:3]

    tri = np.loadtxt(tfile,dtype=np.int32,skiprows=1,comments='#')
    tri[:] = tri[:] - 1
    tri = tri[:,1:4]

    edge = None
    if ( os.path.isfile(efile) ):
        edge = np.loadtxt(efile,skiprows=1,comments='#')
        edge[:,1:3] = edge[:,1:3] - 1
        E = np.array(edge[:,1:4], np.int32)
        bidx = np.where(edge[:,3] != 0)
        bde = E[bidx]
    return pts, tri, bde

def SalvaVTK(arq, p, t, u):
    npts = np.size(p,1)
    nele = np.size(t,1)
    nedg = np.size(e,0)
    
    f = open(arq, 'w')
    f.write("# vtk DataFile Version 2.0\n")
    f.write('2D Unstructured Grid of Linear Triangles\n')
    f.write('ASCII\n')
    f.write('\n')

    f.write('DATASET UNSTRUCTURED_GRID\n')
    f.write('POINTS %d float\n' % (npts))
    for i in range(npts):
        f.write('%f %f %f\n' % (p[0,i],p[1,i],0.0))
    f.write('\n')

    f.write('CELLS %d %d\n' % (nele, 4*nele))
    for i in range(nele):
        f.write('3 %d %d %d\n' % (t[0,i],t[1,i],t[2,i]))
    f.write('\n')

    f.write('CELL_TYPES %d\n' % nele)
    for i in range(nele):
        f.write('5\n')
    f.write('\n')
    
    f.write('POINT_DATA %d\n' % npts)
    f.write('SCALARS scalar_field float\n')
    f.write('LOOKUP_TABLE default\n')
    for i in range(npts):
        f.write('%f\n' % u[i])
    f.write('\n')
    f.close()

# ------------------------------------------------------------------------------    
# Funcoes que dependem do problema - airfoil.1
# ------------------------------------------------------------------------------

def kappa(x,y,bidx):    
    if (x>11.999):
        return 1.0e+10
    else:
        return 0.0

def gD(x,y,bidx):
    return 0.0

def gN(x,y,bidx):   
    if(x<-1.999):
        return 1.0
    else:
        return 0.0
    
def a(x,y):
    return 1.0

def f(x,y):
    return 0.0

# ------------------------------------------------------------------------------

if __name__ == "__main__":


    if(len(sys.argv) != 2):
        print("\n Uso: mef2d_robin arquivo_malha\n")
        sys.exit(1)

    p,t,e = LeMalha(sys.argv[1])
    p = p.transpose()
    t = t.transpose()

    A = MatrizRigidez(p,t,a)
    R = RobinMatrizMassa(p,e,kappa)
    r = RobinVetorCarga(p,e,kappa,gD,gN)
    b = VetorCarga(p,t,f)

    # resolve o sistema
    A = A + R
    r = r + b
    x = np.linalg.solve(A,r)

    print("Solucao")
    print(x)

    # escreve arquivo de dados no formato VTK
    SalvaVTK('saida.vtk',p,t,x)
