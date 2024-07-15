module Plotting

using LinearAlgebra, Plots, LaTeXStrings
using ..RPA.Interactions: interaction
using ..RPA.Response: perform_RPA

export plot_chi, plot_strengthsVsN, plot_QsVsN


function plot_chi(chis::Vector{Matrix{ComplexF64}}; kwargs...)
    xs = kwargs[:path_plot]
    ticks = kwargs[:path_ticks]
    labels = kwargs[:path_labels]

    n = size(chis[begin])[1]
    egs = eigen.(chis)
    egvals = getproperty.(egs, :values)

    p = plot(framestyle = :box, grid=false, title=L"\chi(Q, \Omega=0)")
    for i in 1:n
        plot!(xs, real.(getindex.(egvals, i)), label = "", lw=2.0)
    end
    xticks!(ticks, labels)
    vline!(ticks, color = :black, lw=1.0, ls=:dash, label = "")

    return p
end


function plot_chi(chis::Vector{Matrix{ComplexF64}}, strength::Float64, ks::Vector{Vector{Float64}}; kwargs...)

    interaction_mats = interaction(strength, ks;
        primitives = kwargs[:primitives], subs = kwargs[:subs],
        localDim = kwargs[:localDim], lookup = kwargs[:lookup])

    chi_RPA = perform_RPA.(chis, interaction_mats ; return_matrix = true)

    return plot_chi(chi_RPA; kwargs...)
end


function plot_strengthsVsN(data::Dict, label::String; kwargs...)

    p = plot(framestyle = :box, grid=false, title=label)

    mus = Float64[]
    fillings = Float64[]
    strengths = Float64[]

    for (key, value) in data

        push!(mus, value["mu"])
        push!(fillings, value["filling"])
        push!(strengths, value[label]["critical strength"])
    end

    args = sortperm(mus)

    mus = mus[args]
    fillings = fillings[args]
    strengths = strengths[args]

    plot!(fillings, strengths, label = "", lw=2.0, m=:o)
    xlabel!(L"⟨n⟩")
    ylabel!(L"\mathcal{I}_c")

    return p

end


function plot_QsVsN(data::Dict, label::String; plot_kwargs...)

    p = plot(framestyle = :box, grid=false, title=label, plot_kwargs...)

    mus = Float64[]
    fillings = Float64[]
    Qs = Vector{Float64}[]

    for (key, value) in data

        push!(mus, value["mu"])
        push!(fillings, value["filling"])
        push!(Qs, value[label]["maximum momentum"])
    end

    args = sortperm(mus)

    mus = mus[args]
    fillings = fillings[args]
    Qs = Qs[args]

    plot!(fillings, norm.(Qs), label = "", lw=2.0, m=:o)
    xlabel!(L"⟨n⟩")
    ylabel!(L"|Q_c|")

    return p

end






































end
