using JLD2, TightBindingToolkit, Plots, LaTeXStrings, LinearAlgebra

#####* Define the triangular lattice with 3-site unit cell unit vectors
a1 = 1 / 2 * [ 3.0, sqrt(3)]
a2 = 1 / 2 * [-3.0, sqrt(3)]
l1 = [1.0, 0]
l2 = [0.5, sqrt(3) / 2]
UC = UnitCell([a1, a2], 2, 2)
#####* Define the nearest neighbor distance
const firstNNdistance  = 1.0

#####* Initialize the unit cell with 2 orbitals (spin-up and spin-down) per site
UC = UnitCell( [a1 , a2] , 2 )
AddBasisSite!( UC , 0*l1+0*l2 )
AddBasisSite!( UC , 1*l1+0*l2 )
AddBasisSite!( UC , 0*l1+1*l2 )

#####* The pauli matrices
SpinVec     =   SpinMats(1//2)

#####* Adding a parameter which tracks the nearest neighbor hopping term
const t1  =   -1.0
t1Param   =   Param(t1, 2)
AddIsotropicBonds!(t1Param, UC , firstNNdistance, SpinVec[4] , "t1")

#####* Adding a parameter which tracks the R=1 skyrmion
const J = -10.0
JParam = Param(J, 2)
AddAnisotropicBond!(JParam, UC, 1, 1, [0, 0],  SpinVec[3], 0.0, "skyrmion Hunds")
AddAnisotropicBond!(JParam, UC, 2, 2, [0, 0], -SpinVec[3], 0.0, "skyrmion Hunds")
AddAnisotropicBond!(JParam, UC, 3, 3, [0, 0], -SpinVec[3], 0.0, "skyrmion Hunds")


#####* Adding the parameter to the unit cell
params = [t1Param, JParam]
CreateUnitCell!(UC, params)

Plot_UnitCell!(UC)

const kSize   =   6 * 30 + 3   ##### a Monkhorst grid of size N = 6n+3 covers all the High-symemtry points.
bz            =   BZ([kSize, kSize])    ##### Make the BZ explicitly in 2d
FillBZ!(bz, UC)

path = [bz.HighSymPoints["G"], bz.HighSymPoints["K1"], bz.HighSymPoints["M2"]]

H       =   Hamiltonian(UC, bz)
DiagonalizeHamiltonian!(H)

M = Model(UC, bz, H ; T=0.001, filling=0.5)
SolveModel!(M)

Plot_Band_Structure!(M, path ; labels = [L"\Gamma", L"K", L"M"], plot_legend=false)

#####* Saving the unit cell in a JLD2 file
# file_name = "../../saves/models/triangle_NN.jld2"
# save(file_name, Dict("unit cell" => UC, "parameters" => params))
