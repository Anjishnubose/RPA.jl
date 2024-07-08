using JLD2, TightBindingToolkit, Plots, LaTeXStrings, LinearAlgebra

#####* Define the triangular lattice with 3-site unit cell unit vectors
const l1  = [ 1/2 , sqrt(3)/2 ]
const l2  = [-1/2 , sqrt(3)/2 ]
UC = UnitCell([2*l1, 2*l2], 2, 2)
#####* Define the nearest neighbor distance
const firstNNdistance  = 1.0

#####* Initialize the unit cell with 2 orbitals (spin-up and spin-down) per site
AddBasisSite!( UC , 0*l1+0*l2 )
AddBasisSite!( UC , 1*l1+0*l2 )
AddBasisSite!( UC , 1*l1-1*l2 )
AddBasisSite!( UC , 2*l1-1*l2 )

#####* The pauli matrices
SpinVec     =   SpinMats(1//2)

#####* Adding a parameter which tracks the nearest neighbor hopping term
const V  =   +0.0
VParam   =   Param(V, 2)
AddIsotropicBonds!(VParam, UC , firstNNdistance,
    [1.0;;] ,
    "NN density-density repulsion")

params = [VParam]
values = Dict()
values["NN_repulsive_density-density"] = [1.0]

#####* Saving the unit cell in a JLD2 file
file_name = "../../saves/interactions/triangle_tetra_V.jld2"
save(file_name, Dict("parameters" => params, "values"=>values))
