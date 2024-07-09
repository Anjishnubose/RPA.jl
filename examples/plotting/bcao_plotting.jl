using JLD2, Plots, LaTeXStrings, LinearAlgebra
include("/home/anjishnubose/Research/Repos/RPA.jl/src/RPA.jl")
using .RPA

data = load("/home/anjishnubose/Research/Repos/RPA.jl/saves/data/bcao_Dirac_combined.jld2")
data = data["beta=20.0_mu=0.0"]

thetas = collect(LinRange(0, 1, 201))
labels = ["theta = $(round(theta, digits=3))*2*pi" for theta in thetas]

ks = data["triqs_data"]["contracted"]
ks = Vector{Float64}[eachrow(ks)...]

Js = Float64[]
Qs = Vector{Float64}[]

for label in labels
    push!(Js, data[label]["critical strength"])
    push!(Qs, data[label]["maximum momentum"])
end



function plot_chis(data::Dict, label::String, strength::Float64;
    path_labels::Vector{T} = [L"\Gamma", L"K", L"M", L"\Gamma"]) where {T<:AbstractString}
    chis = data["combined_chis"]
    interaction_mats = data["$(label)_interaction"]
    triqs_data = data["triqs_data"]
    critical_strength = data[label]["critical strength"]

    eigenstates = perform_RPA(chis, strength .* critical_strength .* interaction_mats ; return_matrix=false)
    peak = maxima(eigenstates)

    theta = parse(Float64, label[end-9:end-5])
    ratio = round(tan(theta*2*pi), digits=3)

    p = plot_chi(chis, strength .* critical_strength .* interaction_mats;
        path_plot = triqs_data["path_plot"],
        path_ticks = triqs_data["path_ticks"],
        path_labels = path_labels,
        fontfamily="Computer Modern", guidefontfamily="cm",
        )
    annotate!(8.5, 0.95*peak["maximum eigenvalue"], text(L"J_3/J_1 = %$(ratio)", :black, :center, 12))
    annotate!(8.5, 0.9*peak["maximum eigenvalue"], text(L"J = %$(round(strength * critical_strength, digits=3))", :black, :center, 12))


    savefig("./$(label)_strength=$(strength).png")
    return p

end


########* Plotting peak momentum *#########
function plotQs(data::Dict, Qs::Vector{Vector{Float64}}, thetas::Vector{Float64} ;
        skip_every::Int64 = 1)

    pyplot()

    ks = data["triqs_data"]["contracted"]
    ks = Vector{Float64}[eachrow(ks)...]

    colors = data["triqs_data"]["path_plot"][indexin(Qs, ks)] / data["triqs_data"]["path_plot"][end]

    p = plot(thetas[1:skip_every:end]*pi*2, norm.(Qs[1:skip_every:end]), label = L"|Q_{crit}|",
        framestyle=:box, gridalpha=0.6,
        background_color_legend=RGBA(1.0, 1.0, 1.0, 0.9),
        xticks = ([0.0, 0.25*pi, 0.5*pi, 0.75*pi, pi, 1.25*pi, 1.5*pi, 1.75*pi],
            ["0", L"\frac{\pi}{4}", L"\frac{\pi}{2}", L"\frac{3\pi}{4}", L"\pi", L"\frac{5\pi}{4}", L"\frac{3\pi}{2}", L"\frac{7\pi}{4}"]),
        xlabel = L"\theta=tan^{-1}\left(J_3/J_1\right)",
        marker=:o, lw=2.0, proj=:polar, legend_position=:bottomleft,
        m=cgrad(:darktest, rev=true), zcolor = colors[1:skip_every:end], clims=(0, 1),
        colorbar_ticks = (data["triqs_data"]["path_ticks"]/data["triqs_data"]["path_plot"][end], [L"\Gamma", L"K", L"M", L"\Gamma"]),
        colorbar_title = L"\mathbf{Q}",
        guidefontfamily="cm", fontfamily="cm")

    plot!(thetas*pi*2, repeat([4*pi/3], length(thetas)),
        c=cgrad(:darktest, rev=true)[data["triqs_data"]["path_ticks"][2]/data["triqs_data"]["path_plot"][end]], l=:dash, lw=2,
        label=L"K", linealpha=0.5)

    plot!(thetas*pi*2, repeat([2*pi/sqrt(3)], length(thetas)),
    c=cgrad(:darktest, rev=true)[data["triqs_data"]["path_ticks"][3]/data["triqs_data"]["path_plot"][end]], l=:dash, lw=2,
    label=L"M", linealpha=0.5)

    return p
end

########* Plotting peak strength *#########
function PlotJs(data::Dict, Js::Vector{Float64}, thetas::Vector{Float64} ;
    skip_every::Int64 = 1)

    pyplot()
    ks = data["triqs_data"]["contracted"]
    ks = Vector{Float64}[eachrow(ks)...]

    colors = data["triqs_data"]["path_plot"][indexin(Qs, ks)] / data["triqs_data"]["path_plot"][end]

    p = plot(thetas[1:skip_every:end]*pi*2, (Js[1:skip_every:end]), label = L"|J_{crit}|",
        framestyle=:box, gridalpha=0.6,
        background_color_legend = RGBA(1.0, 1.0, 1.0, 0.4),
        xticks = ([0.0, 0.25*pi, 0.5*pi, 0.75*pi, pi, 1.25*pi, 1.5*pi, 1.75*pi],
            ["0", L"\frac{\pi}{4}", L"\frac{\pi}{2}", L"\frac{3\pi}{4}", L"\pi", L"\frac{5\pi}{4}", L"\frac{3\pi}{2}", L"\frac{7\pi}{4}"]),
        xlabel = L"\theta",
        marker=:o, lw=2.0, proj=:polar, legend_position=:topright,
        m=cgrad(:darktest, rev=true), zcolor = colors[1:skip_every:end], clims=(0, 1),
        colorbar_ticks = (data["triqs_data"]["path_ticks"]/data["triqs_data"]["path_plot"][end], [L"\Gamma", L"K", L"M", L"\Gamma"]),
        colorbar_title = L"\mathbf{Q}",
        guidefontfamily="cm", fontfamily="cm")

    return p
end
