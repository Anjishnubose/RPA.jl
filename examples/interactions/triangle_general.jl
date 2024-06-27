using JLD2, TightBindingToolkit, LinearAlgebra

#####* Define the honeycomb lattice unit vectors
const a1  = [ 1/2 , sqrt(3)/2 ]
const a2  = [-1/2 , sqrt(3)/2 ]
#####* Define the sublattice sites
const b1  = [ 0.0 , 0.0 ]
#####* Define the nearest neighbor and third nearest neighbour distance
const firstNNdistance  = 1.0

#####* Initialize the unit cell with 2 orbitals (spin-up and spin-down) per site
UC = UnitCell( [a1 , a2] , 2 )
AddBasisSite!( UC , b1 )

#####* Adding a parameter which tracks the nearest neighbor Heisenberg interaction in the spin-spin basis
const J1  =   +0.0
J1Param   =   Param(J1, 2)
AddIsotropicBonds!(J1Param, UC , firstNNdistance,
    [0.0 0.0 0.0 0.0 ; 0.0 1.0 0.0 0.0 ; 0.0 0.0 1.0 0.0 ; 0.0 0.0 0.0 1.0] ,
    "J1 Heisenberg")

const U  =   +1.0
UParam   =   Param(U, 2)
AddIsotropicBonds!(UParam, UC , 0.0,
    [1.0 0.0 0.0 0.0 ; 0.0 0.0 0.0 0.0 ; 0.0 0.0 0.0 0.0 ; 0.0 0.0 0.0 -1.0] ,
    "On-site Hubbard")

const V  =   +0.0
VParam   =   Param(V, 2)
AddIsotropicBonds!(VParam, UC , firstNNdistance,
    [1.0 0.0 0.0 0.0 ; 0.0 0.0 0.0 0.0 ; 0.0 0.0 0.0 0.0 ; 0.0 0.0 0.0 0.0] ,
    "NN density-density repulsion")

params = [J1Param, UParam, VParam]

#####* multiple different interactions to run RPA on
values = Dict()
values["NN_AFM_Heisenberg"] = [1.0, 0.0, 0.0]
values["NN_FM_Heisenberg"] = [-1.0, 0.0, 0.0]
values["On-site_repulsive_Hubbard"] = [0.0, 1.0, 0.0]
values["On-site_attractive_Hubbard"] = [0.0, -1.0, 0.0]
values["NN_repulsive_density-density"] = [0.0, 0.0, 1.0]
values["NN_attractive_density-density"] = [0.0, 0.0, -1.0]

#####* Saving the unit cell in a JLD2 file
file_name = "../../saves/interactions/triangle_general.jld2"
save(file_name, Dict("parameters" => params, "values"=>values))
