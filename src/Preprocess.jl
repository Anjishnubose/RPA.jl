module Preprocess

using NPZ, JLD2
using LinearAlgebra

export dress_primitives, dress_reciprocal, dress_chis
export get_reciprocal_ks

#####* Return a vector of the three 3-d real-space primitive vectors from input dictionary
function dress_primitives(data::Dict ; entry::String = "primitives")::Vector{Vector{Float64}}
    primitives = data[entry]
    primitives = Vector{eltype(data[entry])}[eachrow(data[entry])...]
    primitives[1] = [primitives[1]; 0.0]
    primitives[2] = [primitives[2]; 0.0]
    push!(primitives, [0.0, 0.0, 1.0])

    return primitives
end

#####* Returns a vector of the three 3d reciprocal vectors from input dictionary
function dress_reciprocal(data::Dict ; entry::String = "reciprocal")::Vector{Vector{Float64}}
    reciprocal = data[entry]
    reciprocal = Vector{eltype(data[entry])}[eachrow(data[entry])...]
    return reciprocal
end

#####* returns a vector of vectors of all momenta in the BZ in units of the reciprocal lattice vectors.
function get_reciprocal_ks(data::Dict ; entry::String = "bz")::Vector{Vector{Float64}}
    primitives = dress_primitives(data)
    ks = Vector{eltype(data[entry])}[eachrow(data[entry])...]

    k1s = dot.(ks, Ref(primitives[1]))
    k2s = dot.(ks, Ref(primitives[2]))
    k3s = dot.(ks, Ref(primitives[3]))

    Ks = hcat(k1s, k2s, k3s)
    Ks = Vector{eltype(Ks)}[eachrow(Ks)...]

    return Ks
end

#####* returns a vector of (sublattice x spin) susceptibility matrices where the vector corresponds to different momentas.
function dress_chis(data::Dict ; entries::Vector{String} = ["chiXX", "chiYY", "chiZZ"])::Vector{Matrix{ComplexF64}}
    chiXX = data[entries[1]]
    chiYY = data[entries[2]]
    chiZZ = data[entries[3]]

    indX = [1, 4]
    indY = [2, 5]
    indZ = [3, 6]

    chis = Matrix{ComplexF64}[]

    for i in 1:size(chiXX)[1]

        mat = zeros(ComplexF64, 6, 6)

        mat[indX, indX] = chiXX[i, :, :]
        mat[indY, indY] = chiYY[i, :, :]
        mat[indZ, indZ] = chiZZ[i, :, :]

        push!(chis, mat)
    end

    return chis
end




end
