using JLD2, TightBindingToolkit

#####* Define the honeycomb lattice unit vectors
const a1  = [ 1/2 , sqrt(3)/2 ]
const a2  = [-1/2 , sqrt(3)/2 ]
#####* Define the sublattice sites : The two layers
const b1  = [ 0.0 , 0.0 ]
const b2  = [ 0.0 , 1.0/sqrt(3) ]
#####* Define the nearest neighbor distance
const firstNNdistance  = 1.0/sqrt(3)
const secondNNdistance = 1.0

#####* Initialize the unit cell with 2 orbitals (spin-up and spin-down) per site
UC = UnitCell( [a1 , a2] , 2 )
AddBasisSite!( UC , b1 )
AddBasisSite!( UC , b2 )

#####* The pauli matrices
SpinVec     =   SpinMats(1//2)

#####* Adding a parameter which tracks the intra-layer hopping term
const tIntra  =   -1.0
tIntraParam   =   Param(tIntra, 2)
AddIsotropicBonds!(tIntraParam, UC , secondNNdistance, SpinVec[4] , "t-Intra")

#####* Adding the parameter to the unit cell
params = [tIntraParam]
CreateUnitCell!(UC, params)

#####* Saving the unit cell in a JLD2 file
file_name = "../../saves/models/triangle_bilayer.jld2"
save(file_name, Dict("unit cell" => UC, "parameters" => params))
