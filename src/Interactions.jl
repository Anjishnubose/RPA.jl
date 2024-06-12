module Interactions

using LinearAlgebra, TightBindingToolkit

export interaction

#####* returns the interaction matrix at a given momentum k and a lookup table of interaction matrices
function interaction(strength::Float64, k::Vector{Float64},
        primitives::Vector{Vector{Float64}},
        subs::Int64, localDim::Integer,
        lookup::Dict)::Matrix{ComplexF64}

    mat = zeros(ComplexF64, subs*localDim, subs*localDim)

    for (key, value) in lookup
        i, j, offset = key
        δ = sum(offset .* primitives)

        b1  =   localDim * (i - 1) + 1
        b2  =   localDim * (j - 1) + 1

        mat[b1 : b1 + localDim - 1, b2 : b2 + localDim - 1] += value * exp(-1.0im * dot(k, δ))
        mat[b2 : b2 + localDim - 1, b1 : b1 + localDim - 1] += value' * exp(+1.0im * dot(k, δ))
    end

    return strength * mat
end

#####* returns the interaction matrix at a given momentum k and a vector of parameters tracking all the interactions
function interaction(strength::Float64, k::Vector{Float64};
        primitives::Vector{Vector{Float64}},
        subs::Int64, localDim::Integer,
        params::Vector{Param})::Matrix{ComplexF64}

    lookup = Lookup(params)
    return interaction(strength, k,
        primitives, subs, localDim,
        lookup)
end

#####* returns the interaction matrix at a given momentum k and a unit cell of interaction matrices
function interaction(strength::Float64, k::Vector{Float64};
    unitcell::UnitCell)::Matrix{ComplexF64}

    lookup = Lookup(unitcell)
    return interaction(strength, k,
        unitcell.primitives, length(unitcell.basis), unitcell.localDim,
        lookup)
end

#####* returns the interaction matrix at a multiple momenta ks and a lookup table of interaction matrices
function interaction(strength::Float64, ks::Vector{Vector{Float64}},
    primitives::Vector{Vector{Float64}},
    subs::Int64, localDim::Integer,
    lookup::Dict)::Vector{Matrix{ComplexF64}}

    return interaction.(Ref(strength), ks,
        Ref(primitives), Ref(subs), Ref(localDim), Ref(lookup))
    end

#####* returns the interaction matrix at a given momentum k and a vector of parameters tracking all the interactions
function interaction(strength::Float64, ks::Vector{Vector{Float64}};
    primitives::Vector{Vector{Float64}},
    subs::Int64, localDim::Integer,
    params::Vector{Param})::Vector{Matrix{ComplexF64}}

    lookup = Lookup(params)
    return interaction.(Ref(strength), ks,
        Ref(primitives), Ref(subs), Ref(localDim),
        Ref(lookup))
    end

#####* returns the interaction matrix at a given momentum k and a unit cell of interaction matrices
function interaction(strength::Float64, ks::Vector{Vector{Float64}};
unitcell::UnitCell)::Vector{Matrix{ComplexF64}}

    lookup = Lookup(unitcell)
    return interaction.(Ref(strength), ks,
        Ref(unitcell.primitives), Ref(length(unitcell.basis)), Ref(unitcell.localDim),
        Ref(lookup))
    end









end
