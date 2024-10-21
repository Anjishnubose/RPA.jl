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

function anisotropy(D::Float64, E::Float64, theta::Float64)
    return (D*cos(2*theta)-E*sin(2*theta), D*sin(2*theta)+E*cos(2*theta))
end

const λ = 0.0
const D = -0.1
const E = -D
const theta = 2*pi/3

#####* Adding a parameter which tracks the nearest neighbor Heisenberg interaction in the spin-spin basis
const J1  =   +1.0
J1Param   =   Param(J1, 2)
d, e = anisotropy(D, E, 0*theta)
AddAnisotropicBond!(J1Param, UC, 1, 2, [ 0, 0], [1.0+d e 0.0 ; e 1.0-d 0.0 ; 0.0 0.0 λ], firstNNdistance, "J1 XXZ")
d, e = anisotropy(D, E, 1*theta)
AddAnisotropicBond!(J1Param, UC, 1, 2, [-1, 0], [1.0+d e 0.0 ; e 1.0-d 0.0 ; 0.0 0.0 λ], firstNNdistance, "J1 XXZ")
d, e = anisotropy(D, E, 2*theta)
AddAnisotropicBond!(J1Param, UC, 1, 2, [ 0,-1], [1.0+d e 0.0 ; e 1.0-d 0.0 ; 0.0 0.0 λ], firstNNdistance, "J1 XXZ")
# AddIsotropicBonds!(J1Param, UC , firstNNdistance,
#     [1.0 0.0 0.0 ; 0.0 1.0 0.0 ; 0.0 0.0 λ] ,
#     "J1 XXZ")

const J3  =   +1.0
J3Param   =   Param(J3, 2)
AddIsotropicBonds!(J3Param, UC , thirdNNdistance,
    [1.0 0.0 0.0 ; 0.0 1.0 0.0 ; 0.0 0.0 λ] ,
    "J3 XXZ" ; checkOffsetRange = 3)

params = [J1Param, J3Param]

#####* multiple different interactions to run RPA on
const thetas = collect(LinRange(0.0, 1.0, 401))

values = Dict()

for theta in thetas
    values["theta = $(round(theta, digits=3))*2*pi"] = [cos(2*theta*pi), sin(2*theta*pi)]
end

#####* Saving the unit cell in a JLD2 file
file_name = "../../saves/interactions/bcao_J1J3_wAnisotropy.jld2"
save(file_name, Dict("parameters" => params, "values" => values))
