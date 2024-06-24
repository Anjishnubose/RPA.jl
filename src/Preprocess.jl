module Preprocess

using NPZ, JLD2
using LinearAlgebra

export dress_primitives, dress_reciprocal, combine_chis
export get_reciprocal_ks

const labels = Dict(0 => "chi_NN", 1 => "chi_XX", 2 => "chi_YY", 3 => "chi_ZZ", 4 => "chi_NN")

#####* Return a vector of the three 3-d real-space primitive vectors from input dictionary
function dress_primitives(data::Dict ; entry::String = "primitives")::Vector{Vector{Float64}}
    primitives = data[entry]
    primitives = Vector{eltype(data[entry])}[eachrow(data[entry])...]
    # primitives[1] = [primitives[1]; 0.0]
    # primitives[2] = [primitives[2]; 0.0]
    # push!(primitives, [0.0, 0.0, 1.0])

    return primitives
end

#####* Returns a vector of the three 3d reciprocal vectors from input dictionary
function dress_reciprocal(data::Dict ; entry::String = "reciprocal")::Vector{Vector{Float64}}
    reciprocal = data[entry]
    reciprocal = Vector{eltype(data[entry])}[eachrow(data[entry])...]
    return reciprocal
end

#####* returns a vector of vectors of all momenta in the BZ in units of the reciprocal lattice vectors.
function get_reciprocal_ks(data::Dict ; entry::String = "ks")::Vector{Vector{Float64}}
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
function combine_chis(data::Dict ; directions::Vector{Int}, subs::Int64)::Vector{Matrix{ComplexF64}}
    directions = sort(directions)
    data_labels = [labels[i] for i in directions]
    localDim = length(directions)

    chi_combined = [zeros(ComplexF64, subs*localDim, subs*localDim) for _ in 1:size(data[data_labels[begin]])[1]]

    for (ind, label) in enumerate(data_labels)
        chis = data[label]
        inds = [ind + localDim * (i - 1) for i in 1:subs]
        setindex!.(chi_combined, eachslice(chis, dims = 1), Ref(inds), Ref(inds))
    end

    return chi_combined
end




end
