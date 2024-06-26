import numpy as np
#####* triqs utilities
from triqs_tprf.tight_binding import TBLattice
from triqs_tprf.lattice_utils import k_space_path
from triqs_tprf.lattice import lindhard_chi00
from triqs.gf import MeshImFreq, Idx
from triqs_tprf.lattice import lattice_dyson_g0_wk
from triqs_tprf.lattice_utils import imtime_bubble_chi0_wk

#####* Pauli matrices ###########
s1 = np.matrix([[0,1],[1,0]])
s2 = np.matrix([[0,-1j],[1j,0]])
s3 = np.matrix([[1,0],[0,-1]])
s4 = np.eye(2, dtype=np.complex128)
sup = np.matrix([[1,0],[0,0]])
sdo = np.matrix([[0,0],[0,1]])
paulis = [s4, s1, s2, s3, s4]

#####* returns the bare response tensor chi0(k, iwn) for a given model
def bare_chi(beta:float, n_matsubara: int, mu: float, ham):
    
    wmesh = MeshImFreq(beta=beta, S='Fermion', n_max=n_matsubara)
    g0_wk = lattice_dyson_g0_wk(mu=mu, e_k=ham, mesh=wmesh)
    chi00_wk = imtime_bubble_chi0_wk(g0_wk, nw=1)
    return chi00_wk

#####* returns S^a at site i of total sites N where S^a = [rho, Sx, Sy, Sz, rho].
def S(i:int, direction:int, N: int) -> np.matrix:
    mat = np.zeros((2*N, 2*N), dtype=np.complex128)
    mat[2*i:2*i+2, 2*i:2*i+2] = paulis[direction]/2
    
    return mat

#####* contract the full rank-4 susceptibility tensor to return a matrix
def chi_contraction(chi, i:int, j:int, direction:int, N: int):
    Si = S(i, direction, N)
    Sj = S(j, direction, N)
    
    chi_contracted = chi[0, 0, 0, 0].copy()
    chi_contracted.data[:] = np.einsum('wqabcd,ab,cd->wq', chi.data, Si, Sj)[:, :]
    chi_contracted = chi_contracted[Idx(0), :]
    return chi_contracted

#####* interpolate the susceptibility along a high symmetry path in the Brillouin zone
def interpolate_chi(chi_contracted, ks):
    assert( ks.shape[1] == 3 )
    chi_interp = np.zeros(
        [ks.shape[0]] + list(chi_contracted.target_shape), dtype=complex)

    for kidx, (kx, ky, kz) in enumerate(ks):
        chi_interp[kidx] = chi_contracted((kx, ky, kz))

    return chi_interp

#####* interpolate the contracted susceptibility tensor along a high symmetry path in the Brillouin zone
def interpolate_chi_mat(chi, direction:int, N: int, ks):

    chiMats = np.zeros([ks.shape[0], N, N], dtype=complex)
    
    for i in range(N):
        for j in range(N):
            
            chi_contracted = chi_contraction(chi, i, j, direction, N)
            chiMats[:, i, j] += interpolate_chi(chi_contracted, ks)
    
    return chiMats