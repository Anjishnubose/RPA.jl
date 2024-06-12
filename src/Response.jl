module Response

using LinearAlgebra
using ..RPA.Preprocess: dress_chis, get_reciprocal_ks
using ..RPA.Interactions: J_polar, JMats

export perform_RPA, minima, maxima, find_peak

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


function find_peak(theta::Float64, lambda::Float64, data::Dict ;
        steps::Int = 21, lower::Float64 = 0.0, upper::Float64 = 2.0,
        ks_entry::String = "ks")::Dict{String, Any}

    current = Float64[]
    check = nothing
    ##### momentum points.
    reciprocal_ks = get_reciprocal_ks(data ; entry = ks_entry)
    ks = data[ks_entry]
    ks = Vector{eltype(ks)}[eachrow(ks)...]
    ##### bare susceptibility matrices.
    chis = dress_chis(data)
    #####* binary search for the critical |J| at a given theta = tan^{-1}(J3/J1).
    for _ in 1:steps
        push!(current, (upper + lower) / 2)
        ##### determining the interaction matrices.
        J1, J3 = J_polar(current[end], theta)
        Js = JMats(J1, J3, lambda, reciprocal_ks)
        ##### RPA calculation.
        eigenvalues, eigenvectors = perform_RPA(chis, Js)
        check = minima(eigenvalues, eigenvectors)

        if check["minimum eigenvalue"] < 0.0
            upper = current[end]
        else
            lower = current[end]
        end
    end

    if check["minimum eigenvalue"] < 0.0
        J1, J3 = J_polar(lower, theta)
        Js = JMats(J1, J3, lambda, reciprocal_ks)

        eigenvalues, eigenvectors = perform_RPA(chis, Js)
        check = minima(eigenvalues, eigenvectors)
        peak = maxima(eigenvalues, eigenvectors)

        return Dict("J" => lower, check..., maxima...,
                    "minimum reciprocal momentum" => reciprocal_ks[check["minimum index"]],
                    "minimum momentum" => ks[check["minimum index"]],
                    "maximum reciprocal momentum" => reciprocal_ks[peak["maximum index"]],
                    "maximum momentum" => ks[peak["maximum index"]])
    else
        J1, J3 = J_polar(current[end], theta)
        Js = JMats(J1, J3, lambda, reciprocal_ks)

        eigenvalues, eigenvectors = perform_RPA(chis, Js)
        check = minima(eigenvalues, eigenvectors)
        peak = maxima(eigenvalues, eigenvectors)

        return Dict("J" => current[end], check..., maxima...,
                    "minimum reciprocal momentum" => reciprocal_ks[check["minimum index"]],
                    "minimum momentum" => ks[check["minimum index"]],
                    "maximum reciprocal momentum" => reciprocal_ks[peak["maximum index"]],
                    "maximum momentum" => ks[peak["maximum index"]])
    end
end





































































end
