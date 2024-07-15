import numpy as np
import argparse
import yaml
#####* importing other modules
import model as mdl
import bare_response as br

labels = {0 : "chi_NN", 1 : "chi_XX", 2 : "chi_YY", 3 : "chi_ZZ", 4 : "chi_NN"}

if __name__=="__main__":
    
    #####* defining the command line arguments to parse
    parser = argparse.ArgumentParser(
                        prog='ProgramName',
                        description='What the program does',
                        epilog='Text at the bottom of help')
    
    parser.add_argument('input', help='Input file location', type=str, default="")
    args = parser.parse_args()
    #####* loading the input file
    fobj = open(args.input, "r")
    params = yaml.load(fobj, Loader=yaml.CLoader)
    #####* loading the unit cell
    unitcell = np.load(params["unitcell"]["triqs"])
    print("Unit cell loaded")
    
    #####* building the triqs model
    model = mdl.triqs_model(unitcell)
    print("Model built")

    N = int(len(model.orbital_names)/2)
    #####* building the Brillouin zone and a high symmetry path
    ksize = params["k_size"]
    kmesh = model.get_kmesh(n_k=(ksize, ksize, 1))
    ks = np.array([k.value for k in kmesh])
    path_vecs, path_plot, path_ticks = mdl.k_path(model, params["k_points"])
    
    #####* building the hamiltonian
    hamiltonian = mdl.hamiltonian(model, ksize)
    bandwidth = mdl.bandwidth(kmesh, hamiltonian)
    print("Hamiltonian built")
    
    
    #####* fillings vs chemical potential
    beta = params["beta"]
    if "values" in params["mus"]:
        mus = params["mus"]["values"]
    else:
        mus = np.array(np.linspace(*bandwidth, params["mus"]["n"]))
        params["mus"]["values"] = mus.tolist()
        
        with open(args.input, 'w') as file:
            yaml.dump(params, file)
    
    fillings = [mdl.get_filling(mu, beta, hamiltonian, kmesh) for mu in mus]
    
    print("Starting TRIQS calculations...")
    
    for index, mu in enumerate(mus):
        print(f"calculating bare bubble for mu = {mu} => filling = {fillings[index]}...")

        chi00 = br.bare_chi(beta, params["n_matsubara"], mu, hamiltonian)
        
        if params["contract"]=="path":
            ks_contract = path_vecs
        else:
            ks_contract = ks
        
        output = {}
        for direction in params["directions"]:
            chi = br.interpolate_chi_mat(chi00, direction, N, ks_contract)
            output[labels[direction]] = chi

        print("contraction completed")


        fileName = params["output"] + f"_beta={beta}_mu={np.round(mu, 3)}.npz"
        np.savez(fileName, **output,
                    beta = beta, mu = mu, filling=fillings[index], 
                    primitives=model.units, reciprocal = kmesh.bz.units, 
                    ks = ks, path = path_vecs, path_plot = path_plot, path_ticks = path_ticks,
                    contracted = ks_contract, 
                    bandwidth = np.array(bandwidth), bands = np.array([mdl.energies(k, hamiltonian) for k in path_vecs]))
    


