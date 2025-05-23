#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 23 09:53:45 2025

@author: felipe
"""

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt
import copy
from toy_solver import * 
# from numba import njit


def get_ufixed(G, yG, mesh):
    uD = np.zeros((len(mesh.bnd_nodes),2))
    for i, j in enumerate(mesh.bnd_nodes):
        uD[i,:] = G@(mesh.X[j,:] - yG)
        
    return uD


def homogeniseP(mesh, G, u0 = None, component = 'truss'):
    ndofs = 2 * mesh.X.shape[0]
    forces = np.zeros(ndofs)
    yG = np.array([0.5,0.5])
    
    fixed_dofs = np.array([2*mesh.bnd_nodes, 2*mesh.bnd_nodes+1]).T.flatten()
    u_fixed = get_ufixed(G, yG, mesh).flatten()
    
    u0 = np.zeros_like(forces) if type(u0) == type(None) else copy.deepcopy(u0)
    
    u = solve_nonlinear(mesh, forces, fixed_dofs, u_fixed, u0 = u0, component = component)
    
    P = homogeniseP_given_disp(mesh, u, component)
    
    return P, u


def homogeniseP_given_disp(mesh, u, component = 'truss'):
    
    P = np.zeros((2,2))
    stress_list = []
    
    for e in range(mesh.cells.shape[0]):
        n1, n2 = mesh.cells[e]
        dofs = np.array([2*n1, 2*n1+1, 2*n2, 2*n2+1])
        
        XL = mesh.X.flatten()[dofs]
        uL = u[dofs]
        
        Bmat = np.array([[-1,0,1,0], [0,-1,0,1]]) # Delta operation
        a = Bmat@XL
        L0 = np.linalg.norm(a) # underformed length
        a = a/L0 # unitary underformed truss vector 
        Bmat = Bmat/L0 # discrete gradient operator
        
        q = a + Bmat@uL # deformed truss vector (stretch lenght)  
        lmbda = np.linalg.norm(q) # L/L0
        b = q/lmbda # unitary deformed truss vector  
        
        A = mesh.param['A'][e]
        E = mesh.param['E'][e]
        
        V = L0*A
        
        # strain = 0.5*(lmbda**2 - 1)
        # if(component == 'truss'):
        #     stress = E * strain * lmbda
        # elif(component == 'cable'):
        #     stress = E * strain * lmbda if lmbda>1.0 else 0.0

        strain = lmbda - 1
        if(component == 'truss'):
            stress = E * strain 
        elif(component == 'cable'):
            stress = E * strain if lmbda>1.0 else 0.0
            
            
        stress_list.append(A*stress)
        
        P += V*stress*np.outer(b,a)
        
    
    return P


def homogeniseC(mesh, G, component = 'truss'):
    return None 

    # V* a_l D e_k . Bu


# forward finite differences
def homogenise_tang_ffd(mesh, G, tau = 1e-7, component = 'truss'):
    
    P_ref, u_ref = homogeniseP(mesh, G, component = component)
    P_ref = P_ref.flatten()
    Gref = G.flatten()
    n = len(Gref) 
    base_canonic = np.eye(n)
    Atang = np.zeros((n,n))
    
    for j in range(n):
        Gp = (Gref + tau*base_canonic[j,:]).reshape((int(n/2),int(n/2)))
        Pp  = homogeniseP(mesh, Gp, u0 = u_ref, component = component)[0].flatten()
        Atang[:,j] = (Pp - P_ref)/tau 
    
    return Atang

# def get_tangent_pertubation_central_difference(Gmacro, micromodel, tau = 1e-6):
#     n = len(Gmacro)
#     base_canonic = np.eye(n)
#     Atang = np.zeros((n,n))
    
#     for j in range(n):
#         micromodel.restart_initial_guess()
#         micromodel.solve_microproblem(Gmacro + 0.5*tau*base_canonic[j,:])
#         stress_per_p = micromodel.homogenise_stress()
#         micromodel.restart_initial_guess()
#         micromodel.solve_microproblem(Gmacro - 0.5*tau*base_canonic[j,:])
#         stress_per_m = micromodel.homogenise_stress()
        
#         Atang[:,j] = (stress_per_p - stress_per_m)/tau 
    
#     return Atang

