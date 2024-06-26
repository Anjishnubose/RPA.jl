import numpy as np
import yaml
from matplotlib import pyplot as plt
import argparse

labels = {0 : "chi_NN", 1 : "chi_XX", 2 : "chi_YY", 3 : "chi_ZZ", 4 : "chi_NN"}
titles = {0 : r'$\chi_{dd}(\mathbf{Q}, \Omega=0)$', 
            1 : r'$\chi_{xx}(\mathbf{Q}, \Omega=0)$', 
            2 : r'$\chi_{yy}(\mathbf{Q}, \Omega=0)$', 
            3 : r'$\chi_{zz}(\mathbf{Q}, \Omega=0)$', 
            4 : r'$\chi_{dd}(\mathbf{Q}, \Omega=0)$'}

def plot_chi(data: dict, direction: int, path_labels: list, title: str, saveFile: str):
    chi = data[labels[direction]]
    eigs = np.linalg.eigvals(chi)
    
    plt.clf()
    plt.plot(data["path_plot"], np.real(eigs))
    plt.grid()
    plt.xticks(ticks=data["path_ticks"], labels=path_labels)
    plt.title(title)
    plt.savefig(saveFile)



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
    
    #####* fillings vs chemical potential
    beta = params["beta"]
    mus = params["mus"]["values"]
    
    outDirectory = params["output"]
    plotDirectory = params["plots"]
    
    path_labels = [r'${label}$'.format(label = label) for label in params["k_labels"]]
    path_labels.append(r'${label}$'.format(label = params["k_labels"][0]))
    
    for mu in mus:
        output = outDirectory + f"_beta={beta}_mu={np.round(mu, 3)}.npz"
        data = np.load(output)
        
        #####* plotting the bands
        plt.clf()
        plt.plot(data["path_plot"], data["bands"])
        plt.xticks(data["path_ticks"], path_labels)
        plt.ylabel(r'$\epsilon(\mathbf{k})$')
        plt.grid(True)
        plt.savefig(f"{plotDirectory}_bands.png") 
        
        #####* plotting all the bare susceptibilities.
        for direction in params["directions"]:
            plot_title = titles[direction]
            saveFile = f"{plotDirectory}_{labels[direction]}_beta={beta}_mu={np.round(mu, 3)}.png"
            plot_chi(data, direction, path_labels, plot_title, saveFile)
    


