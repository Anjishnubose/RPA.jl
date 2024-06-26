using JLD2, TightBindingToolkit

#####* Define the triangular lattice unit vectors
const a1  = [ 1/2 , sqrt(3)/2 ]
const a2  = [-1/2 , sqrt(3)/2 ]
#####* Define the sublattice sites
const b1  = [ 0.0 , 0.0 ]
#####* Define the nearest neighbor distance
const firstNNdistance  = 1.0

#####* Initialize the unit cell with 2 orbitals (spin-up and spin-down) per site
UC = UnitCell( [a1 , a2] , 2 )
AddBasisSite!( UC , b1 )

#####* The pauli matrices
SpinVec     =   SpinMats(1//2)

#####* Adding a parameter which tracks the nearest neighbor hopping term
const t1  =   -1.0
t1Param   =   Param(t1, 2)
AddIsotropicBonds!(t1Param, UC , firstNNdistance, SpinVec[4] , "t1")

#####* Adding the parameter to the unit cell
params = [t1Param]
CreateUnitCell!(UC, params)

#####* Saving the unit cell in a JLD2 file
file_name = "../../saves/models/triangle_NN.jld2"
save(file_name, Dict("unit cell" => UC, "parameters" => params))
