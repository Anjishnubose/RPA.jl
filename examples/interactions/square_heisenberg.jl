using JLD2, TightBindingToolkit, LinearAlgebra

####* Define the square lattice unit vectors
const a1  = [ 1.0, 0.0 ]
const a2  = [ 0.0, 1.0 ]
#####* Define the sublattice sites
const b1  = [ 0.0 , 0.0 ]
#####* Define the nearest neighbor distance
const firstNNdistance  = 1.0

#####* Initialize the unit cell with 2 orbitals (spin-up and spin-down) per site
UC = UnitCell( [a1 , a2] , 2 )
AddBasisSite!( UC , b1 )

#####* Adding a parameter which tracks the nearest neighbor Heisenberg interaction in the spin-spin basis
const J1  =   +1.0
J1Param   =   Param(J1, 2)
AddIsotropicBonds!(J1Param, UC , firstNNdistance, Matrix{ComplexF64}(I, 3, 3) , "J1")

push!(J1Param.value, -1.0)
params = [J1Param]

#####* Saving the unit cell in a JLD2 file
file_name = "../../saves/interactions/square_heisenberg.jld2"
save(file_name, Dict("parameters" => params))
