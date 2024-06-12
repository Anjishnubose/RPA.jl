import numpy as np
from triqs_tprf.tight_binding import TBLattice
from triqs_tprf.lattice_utils import k_space_path

#####* convert a unitcell dictionary to a triqs model
#####* unitcell dictionary can be made from TightBindingToolkit or from a saved file
def triqs_model(unitcell: dict):
    
    return TBLattice(
        units = unitcell["units"],
        orbital_positions = unitcell["orbital_positions"],
        orbital_names = unitcell["orbital_names"],
        hoppings = unitcell["hoppings"]
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
        
    k_vecs, k_plot, k_ticks = k_space_path(paths, bz=model.bz)
    return k_vecs, k_plot, k_ticks

#####* return bands at a momentum k ##########
def energies(k, ham):
    return np.linalg.eigvalsh(ham(k))

#####* returns band energies at each k-value
def bands(ham, kmesh):
    N = kmesh[:].value.shape[0]
    band = [np.linalg.eigvalsh(ham(kmesh[i].value)) for i in range(N**2)]
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

#####* finding filling vs mu ##################
def filling_vs_mu(beta: float, n: int, ham, kmesh):

    bwidth = bandwidth(kmesh, ham)
    mus = np.linspace(*bwidth, n)
    fillings = np.zeros(len(mus))
    
    print("calculating chemical potential vs filling...")
    for (i, mu) in enumerate(mus):
        fillings[i] = filling(bands(ham, kmesh), beta, mu)
    
    return mus, fillings