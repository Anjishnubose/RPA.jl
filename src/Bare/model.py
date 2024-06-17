import numpy as np
import ast
from triqs_tprf.tight_binding import TBLattice
from triqs_tprf.lattice_utils import k_space_path

#####* convert a unitcell dictionary to a triqs model
#####* unitcell dictionary can be made from TightBindingToolkit or from a saved file
def triqs_model(unitcell: dict):
    
    units = [tuple(unit) for unit in np.transpose(unitcell["units"])]
    positions = [tuple(pos) for pos in np.transpose(unitcell["orbital_positions"])]
    
    subs = len(positions)/2
    names = [f"{i+1}:{spin}" for spin in ["up", "dn"] for i in range(int(subs))]
    hoppings = {tuple(unitcell["hopping offsets"][:, i]) : np.array(unitcell["hopping matrices"][i, :, :]) 
                    for i in range(unitcell["hopping offsets"].shape[1])}
    
    return TBLattice(
        units = units,
        orbital_positions = positions,
        orbital_names = names,
        hoppings = hoppings
    )

#####* return the hamiltonian in the Brillouin zone for the corresponding model
def hamiltonian(model, ksize: int):
    kmesh = model.get_kmesh(n_k=(ksize, ksize, 1))
    ham = model.fourier(kmesh)
    return ham

#####* return a high symmetry path in the Brillouin zone
def k_path(model, k_points: list):
    paths = []
    for i in range(len(k_points)-1):
        paths.append((k_points[i], k_points[i+1]))
    
    paths.append((k_points[-1], k_points[0]))
        
    k_vecs, k_plot, k_ticks = k_space_path(paths, bz=model.bz)
    return k_vecs, k_plot, k_ticks

#####* return bands at a momentum k ##########
def energies(k, ham):
    return np.linalg.eigvalsh(ham(k))

#####* returns band energies at each k-value
def bands(ham, kmesh):
    N = kmesh.dims
    band = [np.linalg.eigvalsh(ham(kmesh[i].value)) for i in range(np.prod(N))]
    band = np.concatenate(band, axis = 0)
    
    return band

#####* return total bandwidth #########
def bandwidth(ks, ham):
    bands = np.array([energies(k, ham) for k in ks])
    return (np.min(bands), np.max(bands))

#####* fermi distribution function
def fermi(e: float, beta: float, mu: float) -> float:
    return 1.0 / (np.exp(beta * (e-mu)) + 1.0)

#####* filling at fixed temperature and chemical potential
def filling(band, beta:float, mu: float)-> float:
    return np.sum(fermi(band, beta, mu))/len(band)

def get_filling(mu: float, beta: float, ham, kmesh) -> float:
    return filling(bands(ham, kmesh), beta, mu)

#####* finding filling vs mu ##################
def filling_vs_mu(beta: float, n: int, ham, kmesh):

    bwidth = bandwidth(kmesh, ham)
    mus = np.linspace(*bwidth, n)
    fillings = np.zeros(len(mus))
    
    print("calculating chemical potential vs filling...")
    for (i, mu) in enumerate(mus):
        fillings[i] = filling(bands(ham, kmesh), beta, mu)
    
    return mus, fillings