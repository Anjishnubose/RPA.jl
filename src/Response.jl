module Response

using LinearAlgebra
using ..RPA.Interactions: interaction

export perform_RPA, minima, maxima, find_instability

#####* Returns the eigenvalues and eigenvectors of the RPA susceptibility matrix at some fixed momentum.
function perform_RPA(chi::Matrix{ComplexF64}, interaction::Matrix{ComplexF64} ;
        return_matrix::Bool = false)
    mat = inv(I + chi * interaction) * chi

    if return_matrix
        return mat
    else
        eigs = eigen(mat)

        values = eigs.values
        vectors = eigs.vectors

        return (values, vectors)
    end
end

#####* Returns the eigenvalues and eigenvectors of the RPA susceptibility matrix at multiple momenta ks.
function perform_RPA(chis::Vector{Matrix{ComplexF64}}, interactions::Vector{Matrix{ComplexF64}} ;
        return_matrix::Bool = false)
    return perform_RPA.(chis, interactions ; return_matrix = return_matrix)
end

#####* returns the eigenvalue and eigenvector corresponding to the minimum eigenvalue over all momenta.
function minima(eigenstates::Vector{Tuple{Vector{ComplexF64}, Matrix{ComplexF64}}})::Dict{String, Any}

    eigenvalues = getindex.(eigenstates, 1)
    eigenvectors = getindex.(eigenstates, 2)

    minEigs = getindex.(eigenvalues, 1)
    value, index = findmin(real.(minEigs))
    return Dict("minimum index" => index,
                "minimum eigenvalue" => value,
                "minimum eigenvector" => eigenvectors[index][:, 1])
end

#####* returns the eigenvalue and eigenvector corresponding to the maximum eigenvalue over all momenta.
function maxima(eigenstates::Vector{Tuple{Vector{ComplexF64}, Matrix{ComplexF64}}})::Dict{String, Any}

    eigenvalues = getindex.(eigenstates, 1)
    eigenvectors = getindex.(eigenstates, 2)

    n = length(eigenvalues[begin])
    minEigs = getindex.(eigenvalues, n)
    value, index = findmax(real.(minEigs))
    return Dict("maximum index" => index,
                "maximum eigenvalue" => value,
                "maximum eigenvector" => eigenvectors[index][:, n])
end


function find_instability(chis::Vector{Matrix{ComplexF64}}, ks::Vector{Vector{Float64}};
        steps::Int = 32, lower::Float64 = 0.0, upper::Float64 = 10.0,
        kwargs...)::Dict{String, Any}

    current = Float64[]
    check = nothing

    #####* binary search for the critical interaction strength |J| at a given unit cell fixing ratios of interactions.
    for _ in 1:steps
        push!(current, (upper + lower) / 2)
        ##### determining the interaction matrices.
        interactions = interaction(current[end], ks ; kwargs...)
        ##### RPA calculation.
        eigenstates = perform_RPA(chis, interactions)

        check = minima(eigenstates)

        if check["minimum eigenvalue"] < 0.0
            upper = current[end]
        else
            lower = current[end]
        end
    end

    primitives = get(kwargs, :primitives, [[0.0, 0.0], [0.0, 0.0]])
    d = length(primitives[begin])

    if check["minimum eigenvalue"] < 0.0
        interactions = interaction(lower, ks ; kwargs...)
        eigenstates = perform_RPA(chis, interactions)

        check = minima(eigenstates)
        peak = maxima(eigenstates)

        k_min = ks[check["minimum index"]]
        k_max = ks[peak["maximum index"]]

        return Dict("critical strength" => lower, check..., peak...,
                    "minimum reciprocal momentum" => dot.(Ref(k_min[1:d]), primitives) ./ (2*pi),
                    "minimum momentum" => k_min,
                    "maximum reciprocal momentum" => dot.(Ref(k_max[1:d]), primitives) ./ (2*pi),
                    "maximum momentum" => k_max)
    else
        interactions = interaction(current[end], ks ; kwargs...)
        eigenstates = perform_RPA(chis, interactions)

        check = minima(eigenstates)
        peak = maxima(eigenstates)

        k_min = ks[check["minimum index"]]
        k_max = ks[peak["maximum index"]]

        return Dict("critical strength" => lower, check..., peak...,
                    "minimum reciprocal momentum" => dot.(Ref(k_min[1:d]), primitives) ./ (2*pi),
                    "minimum momentum" => k_min,
                    "maximum reciprocal momentum" => dot.(Ref(k_max[1:d]), primitives) ./ (2*pi),
                    "maximum momentum" => k_max)
    end
end








































































end
