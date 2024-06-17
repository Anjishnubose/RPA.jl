module Interactions

using LinearAlgebra, TightBindingToolkit

export interaction

#####* returns the interaction matrix at a given momentum k and a lookup table of interaction matrices
function interaction(strength::Float64, k::Vector{Float64};
        primitives::Vector{Vector{Float64}},
        subs::Int64, localDim::Integer,
        lookup::Dict)::Matrix{ComplexF64}

    mat = zeros(ComplexF64, subs*localDim, subs*localDim)

    for (key, value) in lookup
        i, j, offset = key
        δ = sum(offset .* primitives)

        b1  =   localDim * (i - 1) + 1
        b2  =   localDim * (j - 1) + 1

        mat[b1 : b1 + localDim - 1, b2 : b2 + localDim - 1] += value * exp(-1.0im * dot(k[1:length(δ)], δ))
        mat[b2 : b2 + localDim - 1, b1 : b1 + localDim - 1] += value' * exp(1.0im * dot(k[1:length(δ)], δ))
    end

    return strength * mat
end

#####* returns the interaction matrix at a given momentum k and a vector of parameters tracking all the interactions
function interaction(strength::Float64, k::Vector{Float64}, params::Vector{Param};
        primitives::Vector{Vector{Float64}},
        subs::Int64, localDim::Integer)::Matrix{ComplexF64}

    lookup = Lookup(params)
    return interaction(strength, k;
        primitives = primitives, subs = subs, localDim = localDim,
        lookup = lookup)
end

#####* returns the interaction matrix at a given momentum k and a unit cell of interaction matrices
function interaction(strength::Float64, k::Vector{Float64},
    unitcell::UnitCell)::Matrix{ComplexF64}

    lookup = Lookup(unitcell)
    return interaction(strength, k;
        primitives = unitcell.primitives, subs = length(unitcell.basis), localDim = unitcell.localDim,
        lookup = lookup)
end

#####* returns the interaction matrix at a multiple momenta ks and a lookup table of interaction matrices
function interaction(strength::Float64, ks::Vector{Vector{Float64}};
    primitives::Vector{Vector{Float64}},
    subs::Int64, localDim::Integer,
    lookup::Dict)::Vector{Matrix{ComplexF64}}

    return interaction.(Ref(strength), ks;
        primitives = primitives, subs = subs, localDim = localDim,
        lookup = lookup)
    end

#####* returns the interaction matrix at a given momentum k and a vector of parameters tracking all the interactions
function interaction(strength::Float64, ks::Vector{Vector{Float64}}, params::Vector{Param};
    primitives::Vector{Vector{Float64}},
    subs::Int64, localDim::Integer,
    )::Vector{Matrix{ComplexF64}}

    lookup = Lookup(params)
    return interaction.(Ref(strength), ks;
        primitives = primitives, subs = subs, localDim = localDim,
        lookup = lookup)
    end

#####* returns the interaction matrix at a given momentum k and a unit cell of interaction matrices
function interaction(strength::Float64, ks::Vector{Vector{Float64}}, unitcell::UnitCell)::Vector{Matrix{ComplexF64}}

    lookup = Lookup(unitcell)
    return interaction.(Ref(strength), ks;
        primitives = unitcell.primitives, subs = length(unitcell.basis), localDim = unitcell.localDim,
        lookup = lookup)
    end









end
