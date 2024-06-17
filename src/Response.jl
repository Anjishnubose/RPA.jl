module Response

using LinearAlgebra
using ..RPA.Interactions: interaction

export perform_RPA, minima, maxima, find_instability

#####* Returns the eigenvalues and eigenvectors of the RPA susceptibility matrix at some fixed momentum.
function perform_RPA(chi::Matrix{ComplexF64}, interaction::Matrix{ComplexF64})::Tuple{Vector{ComplexF64}, Matrix{ComplexF64}}
    mat = inv(I + chi * interaction) * chi
    eigs = eigen(mat)

    values = eigs.values
    vectors = eigs.vectors

    return (values, vectors)
end

#####* Returns the eigenvalues and eigenvectors of the RPA susceptibility matrix at multiple momenta ks.
function perform_RPA(chis::Vector{Matrix{ComplexF64}}, interactions::Vector{Matrix{ComplexF64}})
    return perform_RPA.(chis, interactions)
end

#####* returns the eigenvalue and eigenvector corresponding to the minimum eigenvalue over all momenta.
function minima(eigenvalues::Vector{Vector{ComplexF64}},
        eigenvectors::Vector{Matrix{ComplexF64}})::Dict{String, Any}

    minEigs = getindex.(eigenvalues, 1)
    value, index = findmin(real.(minEigs))
    return Dict("minimum index" => index,
                "minimum eigenvalue" => value,
                "minimum eigenvector" => eigenvectors[index][:, 1])
end

#####* returns the eigenvalue and eigenvector corresponding to the maximum eigenvalue over all momenta.
function maxima(eigenvalues::Vector{Vector{ComplexF64}},
        eigenvectors::Vector{Matrix{ComplexF64}})::Dict{String, Any}

    n = length(eigenvalues[begin])
    minEigs = getindex.(eigenvalues, n)
    value, index = findmax(real.(minEigs))
    return Dict("maximum index" => index,
                "maximum eigenvalue" => value,
                "maximum eigenvector" => eigenvectors[index][:, n])
end


function find_instability(chis::Vector{Matrix{ComplexF64}}, ks::Vector{Vector{Float64}};
        steps::Int = 21, lower::Float64 = 0.0, upper::Float64 = 2.0,
        kwargs...)::Dict{String, Any}

    current = Float64[]
    check = nothing

    #####* binary search for the critical interaction strength |J| at a given unit cell fixing ratios of interactions.
    for _ in 1:steps
        push!(current, (upper + lower) / 2)
        ##### determining the interaction matrices.
        interactions = interaction(current, ks ; kwargs...)
        ##### RPA calculation.
        eigenvalues, eigenvectors = perform_RPA(chis, interactions)
        check = minima(eigenvalues, eigenvectors)

        if check["minimum eigenvalue"] < 0.0
            upper = current[end]
        else
            lower = current[end]
        end
    end

    if check["minimum eigenvalue"] < 0.0
        interactions = interaction(lower, ks ; kwargs...)
        eigenvalues, eigenvectors = perform_RPA(chis, interactions)

        check = minima(eigenvalues, eigenvectors)
        peak = maxima(eigenvalues, eigenvectors)

        primitives = get(kwargs, "primitives", [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
        d = length(primitives[begin])

        return Dict("critical strength" => lower, check..., maxima...,
                    "minimum reciprocal momentum" => dot.(Ref(ks[check["minimum index"]][1:d]), primitives),
                    "minimum momentum" => ks[check["minimum index"]],
                    "maximum reciprocal momentum" => dot.(Ref(ks[check["maximum index"]][1:d]), primitives),
                    "maximum momentum" => ks[peak["maximum index"]])
    else
        interactions = interaction(current, ks ; kwargs...)
        eigenvalues, eigenvectors = perform_RPA(chis, interactions)

        check = minima(eigenvalues, eigenvectors)
        peak = maxima(eigenvalues, eigenvectors)

        return Dict("critical strength" => current[end], check..., maxima...,
                    "minimum reciprocal momentum" => dot.(Ref(ks[check["minimum index"]][1:d]), primitives),
                    "minimum momentum" => ks[check["minimum index"]],
                    "maximum reciprocal momentum" => dot.(Ref(ks[check["maximum index"]][1:d]), primitives),
                    "maximum momentum" => ks[peak["maximum index"]])
    end
end








































































end
