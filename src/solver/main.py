#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 22 10:32:32 2025

@author: frocha
"""

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt
# from numba import njit

# ----------------------------------------
# Element Stiffness and Internal Force
# ----------------------------------------

# @njit
def truss_element_stiffness_force(X, u, A, E):
    Bmat = np.array([[-1,0,1,0], [0,-1,0,1]]) # Delta operation
    a = Bmat@X
    L0 = np.linalg.norm(a) # underformed length
    a = a/L0 # unitary underformed truss vector 
    Bmat = Bmat/L0 # discrete gradient operator
    
    q = a + Bmat@u # deformed truss vector (stretch lenght)  
    lmbda = np.linalg.norm(q) # L/L0
    b = q/lmbda # unitary deformed truss vector  
    
    # psi = 0.5*E*(lmbda - 1)**2
    strain = lmbda - 1
    stress = E * strain

    V = L0*A  # Volume

    # Internal force vector (2D)
    f_int = V * stress * (Bmat.T @ b)
    D_mat = E * np.outer(b,b) 
    D_geo = stress*(np.eye(2) - np.outer(b,b))/lmbda

    # Tangent stiffness matrix (4x4)
    K = V * Bmat.T@(D_mat + D_geo)@Bmat
    
    return K, f_int

# ----------------------------------------
# Global Assembly
# ----------------------------------------

# @njit
def assemble_global(nnodes, mesh, u):
    ndofs = 2 * nnodes
    data = []
    rows = []
    cols = []
    F_int = np.zeros(ndofs)

    for e in range(elements.shape[0]):
        n1, n2 = mesh.cells[e]
        XL = mesh.X[[n1,n2]].flatten()
        uL = np.array( [ u[2*n1:2*n1+2], u[2*n2:2*n2+2]]).flatten()  
        Ke, fe = truss_element_stiffness_force(XL, uL, 
                                               mesh.param['A'][e], mesh.param['E'][e])

        dofs = np.array([2*n1, 2*n1+1, 2*n2, 2*n2+1])
        F_int[dofs] += fe

        for i in range(4):
            for j in range(4):
                rows.append(dofs[i])
                cols.append(dofs[j])
                data.append(Ke[i, j])

    K_global = sp.coo_matrix((data, (rows, cols)), shape=(ndofs, ndofs)).tocsr()
    return K_global, F_int

# ----------------------------------------
# Apply boundary conditions
# ----------------------------------------

def apply_boundary_conditions(K, F, fixed_dofs, u, u_fixed):
    all_dofs = np.arange(K.shape[0])
    free_dofs = np.setdiff1d(all_dofs, fixed_dofs)

    F_mod = F.copy()
    for i, dof in enumerate(fixed_dofs):
        F_mod -= K[:, dof].toarray().flatten() * u_fixed[i]

    return K[free_dofs][:, free_dofs], F_mod[free_dofs], free_dofs

# ----------------------------------------
# Newton-Raphson Solver
# ----------------------------------------

def solve_nonlinear_truss(mesh, forces, fixed_dofs, u_fixed, tol=1e-6, max_iter=30):
    nnodes = coords.shape[0]
    ndofs = 2 * nnodes
    u = np.zeros(ndofs)

    for iter in range(max_iter):
        K, F_int = assemble_global(nnodes, mesh, u)
        R = forces - F_int

        K_mod, R_mod, free_dofs = apply_boundary_conditions(K, R, fixed_dofs, u, u_fixed)
    
        du = np.zeros_like(u)
        du[free_dofs] = spla.spsolve(K_mod, R_mod)

        u += du

        norm_res = np.linalg.norm(R_mod)
        print(f"Iter {iter:2d}: Residual = {norm_res:.3e}")
        if norm_res < tol:
            break
    return u

# plotting
def plot_truss(coords, elements, u, scale=1.0, show_nodes=True):
    """
    Plot undeformed and deformed truss structure.

    Parameters:
        coords (ndarray): Original coordinates (n_nodes x 2)
        elements (ndarray): Element connectivity (n_elements x 2)
        u (ndarray): Global displacement vector (2*n_nodes)
        scale (float): Scale factor for displacements
        show_nodes (bool): Show node indices
    """
    n_nodes = coords.shape[0]
    u_nodes = u.reshape((n_nodes, 2))
    coords_def = coords + scale * u_nodes

    plt.figure(figsize=(8, 6))
    for e in elements:
        n1, n2 = e
        x_orig = coords[[n1, n2]]
        x_def = coords_def[[n1, n2]]

        # Undeformed (dashed gray)
        plt.plot(x_orig[:, 0], x_orig[:, 1], 'k--', lw=1, alpha=0.5)

        # Deformed (solid red)
        plt.plot(x_def[:, 0], x_def[:, 1], 'r-', lw=2)

    if show_nodes:
        for i, (x, y) in enumerate(coords):
            plt.text(x, y, f'{i}', color='blue', fontsize=10)

    plt.axis('equal')
    plt.grid(True)
    plt.title("Truss: Undeformed (black dashed) and Deformed (red)")
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.show()

class Mesh:
    def __init__(self, X, cells, param=None):
        self.X = X
        self.cells = cells
        self.param = param
        self.ndim = self.X.shape[1]
        

if __name__ == "__main__":
    # Node coordinates
    coords = np.array([
        [0.0, 0.0],
        [1.0, 0.0],
        [0.5, 1.0],
    ])
    # Elements (0-based indexing)
    elements = np.array([
        [0, 2],
        [1, 2],
        [0, 1],
    ])
    
    A = np.array([1e-4, 1e-4, 1e-4])
    E = np.array([210e9, 210e9, 210e9])

    
    mesh = Mesh(coords, elements, param = {'E': E, 'A': A})

    ndofs = 2 * coords.shape[0]
    forces = np.zeros(ndofs)
    forces[5] = -1e6  # Load at node 2, y-direction

    fixed_dofs = np.array([0, 1, 2, 3])  # Node 0 and 1 fixed
    u_fixed = np.zeros_like(fixed_dofs, dtype=float)

    u = solve_nonlinear_truss(mesh, forces, fixed_dofs, u_fixed)

    print("\nDisplacements:")
    print(u.reshape(-1, 2))
    
    plot_truss(coords, elements, u, scale=10.0)
