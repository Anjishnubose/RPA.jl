using JLD2, TightBindingToolkit, LinearAlgebra

#####* Define the triangular lattice unit vectors
const a1  = [ 1/2 , sqrt(3)/2 ]
const a2  = [-1/2 , sqrt(3)/2 ]
#####* Define the sublattice sites
const b1  = [ 0.0 , 0.0 ]
#####* Define the nearest neighbor distance


function get_distances(primitives::Vector{Vector{Float64}}, check_range::Int64)

    distances = [norm(sum([i, j] .* primitives)) for i in -check_range:check_range, j in -check_range:check_range]
    distances = reshape(distances, length(distances))
    sort!(distances)

    bonds = Int64[1]
    unique_distances = [distances[1]]
    k=1
    for (d, distance) in enumerate(distances[2:end])
        if abs(distance-distances[d])<1e-6
            bonds[k] +=1
        else
            push!(bonds, 1)
            push!(unique_distances, distance)
            k+=1
        end
    end

    return (unique_distances, bonds)
end

const firstNNdistance  = 1.0
distances, counts = get_distances([a1, a2], 10)

#####* Initialize the unit cell with 2 orbitals (spin-up and spin-down) per site
UC = UnitCell( [a1 , a2] , 2 )
AddBasisSite!( UC , b1 )

#####* The pauli matrices
SpinVec     =   SpinMats(1//2)

#####* Adding a parameter which tracks the nearest neighbor hopping term
HoppingParams = Param{2, Float64}[]
for (d, distance) in enumerate(distances[2:end])

    if counts[d+1] % 3 ==0
        param   =   Param(-1.0, 2)
        AddIsotropicBonds!(param, UC , distance, SpinVec[4] / distance , "t_$(d)" ; checkOffsetRange = 10)
        push!(HoppingParams, param)
    end

end
# const t1  =   -1.0
# t1Param   =   Param(-1.0, 2)
# AddIsotropicBonds!(t1Param, UC , firstNNdistance, SpinVec[4] / firstNNdistance , "t1")


const kSize   =   6 * 30 + 3   ##### a Monkhorst grid of size N = 6n+3 covers all the High-symemtry points.
bz            =   BZ([kSize, kSize])    ##### Make the BZ explicitly in 2d
FillBZ!(bz, UC)




#####* Adding the parameter to the unit cell
# params = [t1Param]
CreateUnitCell!(UC, HoppingParams)

#####* Saving the unit cell in a JLD2 file
# file_name = "../../saves/models/triangle_NN.jld2"
# save(file_name, Dict("unit cell" => UC, "parameters" => params))

path = [bz.HighSymPoints["G"], bz.HighSymPoints["K1"], bz.HighSymPoints["M2"]]

H       =   Hamiltonian(UC, bz)
DiagonalizeHamiltonian!(H)

M = Model(UC, bz, H ; T=0.001, filling=0.5)
SolveModel!(M)
