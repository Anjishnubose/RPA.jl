using JLD2, TightBindingToolkit, LinearAlgebra

#####* Define the honeycomb lattice unit vectors
const a1  = [ 1/2 , sqrt(3)/2 ]
const a2  = [-1/2 , sqrt(3)/2 ]
#####* Define the sublattice sites
const b1  = [ 0.0 , 0.0 ]
const b2  = [ 0.0 , 1.0/sqrt(3) ]
#####* Define the nearest neighbor and third nearest neighbour distance
const firstNNdistance  = 1.0/sqrt(3)
const secondNNdistance = 1.0

#####* Initialize the unit cell with 2 orbitals (spin-up and spin-down) per site
UC = UnitCell( [a1 , a2] , 2 )
AddBasisSite!( UC , b1 )
AddBasisSite!( UC , b2 )

#####* Adding a parameter which tracks the inter-layer density-density interaction
const U  =   +1.0
UParam   =   Param(U, 2)
AddAnisotropicBond!(UParam, UC, 1, 2, [0, 0], [1.0;;], firstNNdistance, "inter-layer density-density")

#####* Adding a parameter which tracks the intra-layer density-density interaction
const V  =   +1.0
VParam   =   Param(V, 2)
AddIsotropicBonds!(VParam, UC , secondNNdistance,
    [1.0;;] ,
    "intra-layer density-density")

params = [VParam, UParam]

#####* multiple different interactions to run RPA on
values = Dict()
const thetas = collect(LinRange(0, 1, 51))

for theta in thetas
    values["theta = $(round(theta/4, digits=3))*2*pi"] = [cos(theta*pi/2), sin(theta*pi/2)]
end




#####* Saving the unit cell in a JLD2 file
file_name = "../../saves/interactions/triangle_bilayer_UV.jld2"
save(file_name, Dict("parameters" => params, "values" => values))
