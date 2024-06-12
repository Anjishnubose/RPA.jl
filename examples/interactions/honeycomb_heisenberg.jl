using JLD2, TightBindingToolkit, LinearAlgebra

#####* Define the honeycomb lattice unit vectors
const a1  = [ 1/2 , sqrt(3)/2 ]
const a2  = [-1/2 , sqrt(3)/2 ]
#####* Define the sublattice sites
const b1  = [ 0.0 , 0.0 ]
const b2  = [ 0.0 , 1.0/sqrt(3) ]
#####* Define the nearest neighbor and third nearest neighbour distance
const firstNNdistance  = 1.0/sqrt(3)
const thirdNNdistance  = 2/sqrt(3)

#####* Initialize the unit cell with 2 orbitals (spin-up and spin-down) per site
UC = UnitCell( [a1 , a2] , 2 )
AddBasisSite!( UC , b1 )
AddBasisSite!( UC , b2 )

#####* Adding a parameter which tracks the nearest neighbor Heisenberg interaction in the spin-spin basis
const J1  =   -1.0
J1Param   =   Param(J1, 2)
AddIsotropicBonds!(J1Param, UC , firstNNdistance, Matrix{ComplexF64}(I, 3, 3) , "J1")

#####* Adding a parameter which tracks the third nearest neighbor Heisenberg interaction in the spin-spin basis
const J3  =   +0.25
J3Param   =   Param(J3, 2)
AddIsotropicBonds!(J3Param, UC , thirdNNdistance, Matrix{ComplexF64}(I, 3, 3) , "J3")

#####* Adding the parameter to the unit cell
params = [J1Param, J3Param]
CreateUnitCell!(UC, params)

#####* Saving the unit cell in a JLD2 file
file_name = "../../saves/interactions/honeycomb_heisenberg.jld2"
save(file_name, Dict("unit cell" => UC, "parameters" => params))
