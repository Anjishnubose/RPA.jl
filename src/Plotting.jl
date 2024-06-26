module Plotting

using LinearAlgebra, Plots, LaTeXStrings
using ..RPA.Interactions: interaction
using ..RPA.Response: perform_RPA

export plot_chi


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












































end
