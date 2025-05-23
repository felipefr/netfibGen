#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 21 10:44:13 2025

@author: frocha
"""

import numpy as np
import os, sys

sys.path.append('../../src/')
#~ import matrixAndFibresGeneration as mfg
#~ import newTheoryFibresLib as ntf 
#from matrixAndFibresGeneration import * 
#from newTheoryFibresLib import * 


import copy
from fibresLib import *
from toy_solver import *
from toy_homogeniser import *


def get_network():
    np.random.seed(5)
    
    net = Network()
    
    maxit = 99999
    maxnochange = 100
    dJmin = 0.01
    Jtol = 1.0e-10
    ampP0 = 0.01
    alphaP = 0.9
    gammaP = 0.95 # damping max coordinate, for min and max. 0 means no damping (fixed window) and 1 means maximum damping (unit window).
    ampA0 = 0.01
    maxA = 1.1
    minA = 0.9
    alphaA = 0.9
    pertPoint = 0.01 # 0.0 for regular
    alphaPert = 0.5
    restartPert = 4
    ksmooth = 200
    omegaSmooth = 0.5
    timesSmooth = 4
    Jtol2 = 1.0e-15 
    
    setParamOpt = [maxit, maxnochange, dJmin, Jtol, ampP0, alphaP, gammaP, ampA0, alphaA, maxA, minA, pertPoint, alphaPert, restartPert, omegaSmooth, timesSmooth, Jtol2]
    
    pertPoint = 0.0 # 0.0 for regular
    
    net.asymFac = -2 # nx - ny
    net.nFibPrevision = 90
    
    net.createNetwork()
    net.removeVertHoriFibers()
    
    net.setFlagsAndConnectivity()
    net.set_lfa(2,[1.00,0.0])
    net.setAf(2,[0.1,0.0])
    
    
    Preg1 = copy.deepcopy(net.P)
    net.addPertubation(pertPoint)
    net.correctPoints(0.99,0.01)
    
    net.set_af_Lf_Vf()
    net.setNormalAndAbar()
    
    #~ print net.Abar
    #~ print net.normal
    
    AfOld = copy.deepcopy(net.Af)
    
    colour = 'black'
    writeFigNetwork(net,c= colour, figNum = 2 , filename = 'networkNotOptimised.pdf')
    
    #~ net.optimize(setParamOpt,functionalNBCfib)
    #~ net.optimize(setParamOpt,functionalNBCnormal)
    net.optimize(setParamOpt,functionalNBCBoth)
    
    #~ writeFigNetwork(net,figNum = 2, filename = 'networkOptimised.png')
    
    #Ndof = 6
    #Nsubsteps = 6
    #Nparam = 29
    # net.writeNetwork(Ndof,Nsubsteps,Nparam,fignum = 1,opInifile = 'Piola', opIncludeTri = 0, addAuxNodes = 2, addAuxNodeOnFibres = 1)
    
    return net

# =========== using the toy solver ====================

net = get_network()
# mesh = Mesh(net.P, net.ElemFib, param = {'E': len(net.ElemFib)*[100.0] , 'A' : net.Af + 0.1*np.random.rand(len(net.Af))})

mesh = Mesh(net.P, net.ElemFib, param = {'E': len(net.ElemFib)*[100.0] , 
                                         'A' : net.Af + 0.1*np.random.rand(len(net.Af)),
                                         # 'A' : net.Af,
                                         'eta': len(net.ElemFib)*[0.0001]})

mesh.mark_boundary_nodes()


tlist = np.linspace(0.0,1.0,10)
Plist1 = []
Plist2 = []
Clist1 = []
Clist2 = []

Gmax = np.array([[0.15,0.15],[0.1,-0.1]])


for i, t in enumerate(tlist):
    G = t*Gmax
    P, uold = homogeniseP(mesh, G, component = 'truss')
    Plist1.append(P)
    Clist1.append(homogeniseC(mesh, G, component = 'truss'))

    

for i, t in enumerate(tlist):
    G = t*Gmax
    u_fixed = get_ufixed(G).flatten()
    
    u0 = np.zeros_like(forces) if i==0 else u2
    
    u2 = solve_nonlinear(mesh, forces, fixed_dofs, u_fixed, u0 = u0, component = 'cable')
    Plist2.append(homogeniseP(mesh, u2, component = 'cable'))
    Clist2.append(homogeniseC(mesh, u2, component = 'cable'))


Paux, force_list = homogeniseP(mesh, u2, component = 'cable')

plot_truss(mesh, u1, scale=1.0)
plot_truss(mesh, u2, scale=1.0)

# plot_truss(mesh, u2, scale=1.0)

Plist1 = np.array(Plist1)
Plist2 = np.array(Plist2)

plt.title('homogenised stress P11')
# plt.plot(tlist, Plist1[:,1,1], '-o', label = 'truss')
# plt.plot(tlist, Plist2[:,1,1], '-o', label = 'cables')
plt.plot(tlist, Plist1[:,1,1], '-o', label = 'truss')
plt.plot(tlist, Plist2[:,1,1], '-o', label = 'cables')
plt.legend()
plt.grid()

# print(np.linalg.norm(u1-u2))







