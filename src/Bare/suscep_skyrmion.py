from triqs_tprf.tight_binding import TBLattice
import numpy as np
import matplotlib.pyplot as plt
from triqs_tprf.lattice_utils import k_space_path
import inspect
import string
from triqs_tprf.lattice import lindhard_chi00
from triqs.gf import MeshImFreq, Idx
from triqs_tprf.lattice import lattice_dyson_g0_wk
from triqs_tprf.lattice_utils import imtime_bubble_chi0_wk

import argparse
parser = argparse.ArgumentParser(
                    prog='ProgramName',
                    description='What the program does',
                    epilog='Text at the bottom of help')
parser.add_argument('B', help='Hunds coupling', type=float, default=0.0)
args = parser.parse_args()

loc = "./"
###### Pauli matrices ###########
s1 = np.matrix([[0,1],[1,0]])
s2 = np.matrix([[0,-1j],[1j,0]])
s3 = np.matrix([[1,0],[0,-1]])
s4 = np.eye(2, dtype=np.complex128)
sup = np.matrix([[1,0],[0,0]])
sdo = np.matrix([[0,0],[0,1]])
paulis = [s4, s1, s2, s3]

##### *return the Hunds coupling matrix
def hunds(B: list) -> np.matrix:
    return (B[0]*s1 + B[1]*s2 + B[2]*s3)/2

##### *Triangular Lattice primitives
l1 = np.array([1.0, 0.0, 0.0])
l2 = np.array([-0.5, np.sqrt(3)/2, 0.0])

##### *returns the positions of the atoms in the triangular lattice for a n1 x n2 supercell
def positions(n1: int, n2:int, orbitals:int) -> list:
    pos = []
    names = []
    orbital_names = list(string.ascii_lowercase[0:orbitals])
    for j in range(n2):
        for i in range(n1):
            for orbital in range(orbitals):
                pos.append((i*l1[0] + j*l2[0], i*l1[1] + j*l2[1], 0.0))
                names.append(str(n2*i+j+1) + orbital_names[orbital])
    return pos, names
    
######### *returns the correct hunds coupled matrix for a skyrmion of size=2 ##############
def skyrmion():
    B_center = [0.0, 0.0, 1.0]
    B_edge = [0.0, 0.0, -1.0]
    B_mid1 = [-1.0, 0.0, 0.0]
    B_mid2 = [-0.5, np.sqrt(3)/2, 0.0]
    B_mid3 = [0.5, np.sqrt(3)/2, 0.0]
    B_mid4 = [1.0, 0.0, 0.0]
    B_mid5 = [0.5, -np.sqrt(3)/2, 0.0]
    B_mid6 = [-0.5, -np.sqrt(3)/2, 0.0]
    
    Bs = {"center" : B_center, "edge" : B_edge, 
                "mid1" : B_mid1, "mid2" : B_mid2, "mid3" : B_mid3, 
                "mid4" : B_mid4, "mid5" : B_mid5, "mid6" : B_mid6}
    
    positions = {12 : "edge", 
                        1 : "center", 10 : "mid6", 5 : "edge", 11 : "mid5", 8 : "mid3",
                        4 : "edge", 9 : "edge", 3 : "edge", 6: "mid1", 7 : "mid2", 2 : "mid4"}
    
    mat = np.zeros((24, 24), dtype=np.complex128)
    
    for i in range(12):
        mat[2*i:2*i+2, 2*i:2*i+2] += hunds(Bs[positions[i+1]])
        
    return mat
    
#######* returns the nearest neighbour hopping matrix ##########
def hopping(orbital_space: np.matrix, neighbours: list) -> np.matrix:
    mat = np.zeros((24, 24), dtype=np.complex128)
    for neighbour in neighbours:
        i,j = neighbour
        mat[2*(i-1):2*(i-1)+2, 2*(j-1):2*(j-1)+2] += orbital_space
    return mat

#####* Nearest neighbours within unit cell : 21
NN_00 = [(1, 2), (1, 7), (1, 8), (2, 3), (2, 8), (2, 9), (3, 4), (3, 9), (3, 10), 
        (4, 5), (4, 10), (4, 11), (5, 6), (5, 11), (5, 12), 
        (6, 12), 
        (7, 8), (8, 9), (9, 10), (10, 11), (11, 12)] 
#####* between unit cell at (1, 0) : 7
NN_10 = [(7, 3), (7, 4), (8, 4), (8, 5), (9, 5), (9, 6), (10, 6)]
#####* between unit cell at (0, 1) : 5
NN_01 = [(10, 1), (11, 1), (11, 2), (12, 2), (12, 3)] 
#####* between unit cell at (1, 1) : 1
NN_1m1 = [(6, 1), (6, 7), (12, 7)]
#####* all bonds
bonds = {(0, 0) : NN_00, (0, 1) : NN_01, (1, 0) : NN_10, (-1, 1) : NN_1m1}

#####* Fixed hopping parameters
t1 = 1.0
B = args.B * t1    ##### *Hund's coupling
L = (6, 2)    #####* supercell size
N = L[0] * L[1]

pos, names = positions(*L, 2)
##########* tight binding model #################
model = TBLattice(
    units = [(-3.0, np.sqrt(3), 0.0), [3.0, np.sqrt(3), 0.0]],##### Primitive vectors
    orbital_positions = pos,
    orbital_names = names,
    hoppings = {
        ( 0, 0): -t1 * hopping(np.eye(2), bonds[(0, 0)]) - t1 * np.conj(np.transpose(hopping(np.eye(2), bonds[(0, 0)]))) - B * skyrmion(),##### within the unit cell
        (+1, 0): -t1 * hopping(np.eye(2), bonds[(1, 0)]),
        (-1, 0): -t1 * np.conj(np.transpose(hopping(np.eye(2), bonds[(1, 0)]))),
        ( 0,+1): -t1 * hopping(np.eye(2), bonds[(0, 1)]),
        ( 0,-1): -t1 * np.conj(np.transpose(hopping(np.eye(2), bonds[(0, 1)]))),
        (-1, 1): -t1 * hopping(np.eye(2), bonds[(-1, 1)]),
        ( 1,-1): -t1 * np.conj(np.transpose(hopping(np.eye(2), bonds[(-1, 1)]))),

        },
    )

#####* BZ ######
ksize = 32
kmesh = model.get_kmesh(n_k=(ksize, ksize, 1))
e_k = model.fourier(kmesh)
#####* Path in the BZ
G = [0.0, 0.0, 0.0]
K1 = [1/3, 1/3, 0.0]
K2 = [2/3, 2/3, 0.0]
M1 = [0.5, 0.0, 0.0]
M2 = [0.0, 0.5, 0.0]
M3 = [0.5, 0.5, 0.0]

paths = [(G, K1), (K1, M2), (M2, G)]
k_vecs, k_plot, k_ticks = k_space_path(paths, bz=model.bz)

#####* return bands at a momentum k ##########
def energies(k, e_k):
    return np.linalg.eigvalsh(e_k(k))

#####* return total bandwidth #########
def bandwidth(k_vecs, e_k):
    bands = np.array([energies(k, e_k) for k in k_vecs])
    return (np.min(bands), np.max(bands))

#####* Band structure ##############################
bands_plot = plt.plot(k_plot, [energies(k, e_k) for k in k_vecs])
plt.xticks(k_ticks, [r'$\Gamma$',r'$K_1$',r'$M_2$',r'$\Gamma$'])
plt.ylabel(r'$\epsilon(\mathbf{k})$')
plt.grid(True)
plt.savefig(loc+f"band_structure_t1={t1}_B={B}.png")
plt.close()

quit()

##### fermi distribution function
def fermi(e, beta, mu):
    return 1.0 / (np.exp(beta * (e-mu)) + 1.0)
##### returns band energies at each k-value
def bands(e_k, kmesh, N):
    band = [np.linalg.eigvalsh(e_k(kmesh[i].value)) for i in range(N**2)]
    band = np.concatenate(band, axis = 0)
    
    return band

#####* filling at fixed temperature and chemical potential
def filling(band, beta, mu):
    return np.sum(fermi(band, beta, mu))/len(band)

#####* finding filling vs mu ##################
beta = 20.0 ##### Temperature
bandwidth = bandwidth(k_vecs, e_k)
mus = np.linspace(*bandwidth, 41)
fillings = np.zeros(len(mus))
print("calculating chemical potential vs filling...")
for (i, mu) in enumerate(mus):
    fillings[i] = filling(bands(e_k, kmesh, ksize), beta, mu)
    
mu_vs_filling = plt.plot(fillings, mus)
plt.ylabel(r'$\mu$')
plt.xlabel(r'$n$')
plt.vlines(0.5, *bandwidth, ls='--', color='orange')
plt.savefig(loc + f"mu_vs_filling_t1={t1}_B={B}.png")
plt.close()

############################################* BARE BUBBLE #######################################
#####* returns S^a at site i of total sites N where S^a = [rho, Sx, Sy, Sz].
def S(i:int, N:int, direction:int) -> np.matrix:
    mat = np.zeros((2*N, 2*N), dtype=np.complex128)
    mat[2*i:2*i+2, 2*i:2*i+2] = paulis[direction]
    return mat

#####* contract the full rank-4 susceptibility tensor to return a matrix
def chi_contraction(chi, i:int, j:int, N:int, direction:int):
    Si = S(i, N, direction)
    Sj = S(j, N, direction)
    
    chi_density = chi[0, 0, 0, 0].copy()
    chi_density.data[:] = np.einsum('wqabcd,ab,cd->wq', chi.data, Si, Sj)[:, :]
    chi_density = chi_density[Idx(0), :]
    return chi_density

def interpolate_chi(chi, k_vecs):
    assert( k_vecs.shape[1] == 3 )
    chi_interp = np.zeros(
        [k_vecs.shape[0]] + list(chi.target_shape), dtype=complex)

    for kidx, (kx, ky, kz) in enumerate(k_vecs):
        chi_interp[kidx] = chi((kx, ky, kz))

    return chi_interp

def interp_mat(chi, N:int, direction:int):

    chiMats = np.zeros([k_vecs.shape[0], N, N], dtype=complex)
    
    for i in range(N):
        for j in range(i+1):
            
            chi_contracted = chi_contraction(chi, i, j, N, direction)
            chiMats[:, i, j] += interpolate_chi(chi_contracted, k_vecs)

            if i != j:
                chiMats[:, j, i] += chiMats[:, i, j].conj() 
    
    return chiMats

#####* Mesh for the bubble ##########
wmesh = MeshImFreq(beta=beta, S='Fermion', n_max=60)
ks = np.array([k.value for k in kmesh])

from triqs_tprf.lattice_utils import imtime_bubble_chi0_wk

for index, mu in enumerate(mus):
    print(f"calculating bare bubble for mu = {mu} => filling = {fillings[index]}...")

    g0_wk = lattice_dyson_g0_wk(mu=mu, e_k=e_k, mesh=wmesh)
    chi00_wk = imtime_bubble_chi0_wk(g0_wk, nw=1)
    print("TRIQS calculation done!")
    chi_density = interp_mat(chi00_wk, N, 0)

    # chiDD = np.zeros((k_vecs.shape[0], N, N), dtype = np.complex128)

    # for i in range(N):
    #     for j in range(i+1):
    #         chiDD[:, i, j] = chi_contraction(chi00_wk, i, j, N, 0).data
            
    #         if i != j:
    #             chiDD[:, j, i] = np.conj(chiDD[:, i, j])
    print("contraction completed")

    fileName = f"t1={t1}_B={B}_beta={beta}_mu={mu}_suscep.npz"
    np.savez(loc+fileName, chiDD = chi_density, 
                ks=k_vecs, beta = beta, mu = mu, filling=fillings[index],
                reciprocal = kmesh.bz.units, primitives=model.units)