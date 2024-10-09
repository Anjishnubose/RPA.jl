using JLD2, Plots, LaTeXStrings, LinearAlgebra
include("/home/anjishnubose/Research/Repos/RPA.jl/src/RPA.jl")
using .RPA



function plot_chis(data::Dict, label::String, strength::Float64;
    path_labels::Vector{T} = [L"\Gamma", L"K", L"M", L"\Gamma"]) where {T<:AbstractString}

    chis = data["combined_chis"]
    interaction_mats = data["$(label)_interaction"]
    triqs_data = data["triqs_data"]
    critical_strength = data[label]["critical strength"]

    eigenstates = RPA.perform_RPA(chis, strength .* critical_strength .* interaction_mats ; return_matrix=false)
    peak = RPA.maxima(eigenstates)

    label_split = split(label, "=")
    label_split = split(label_split[end], "*")
    theta = parse(Float64, label_split[begin])
    ratio = round(tan(theta*2*pi), digits=3)

    p = RPA.plot_chi(chis, strength .* critical_strength .* interaction_mats;
        path_plot = triqs_data["path_plot"],
        path_ticks = triqs_data["path_ticks"],
        path_labels = path_labels,
        fontfamily="Computer Modern", guidefontfamily="cm",
        )
    annotate!(8.5, 0.95*peak["maximum eigenvalue"], text(L"J_3/J_1 = %$(ratio)", :black, :center, 12))
    annotate!(8.5, 0.9*peak["maximum eigenvalue"], text(L"J = %$(round(strength * critical_strength, digits=3))", :black, :center, 12))


    # savefig("./$(label)_strength=$(strength).png")
    return p

end


########* Plotting peak momentum *#########
function PlotQs(data::Dict, Qs::Vector{Vector{Float64}}, thetas::Vector{Float64} ;
        skip_every::Int64 = 1, title::AbstractString = "")

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
        markersize = 5, markerstrokealpha=0.25,
        m=cgrad(:darktest, rev=true), zcolor = colors[1:skip_every:end], clims=(0, 1),
        colorbar_ticks = (data["triqs_data"]["path_ticks"]/data["triqs_data"]["path_plot"][end], [L"\Gamma", L"K", L"M", L"\Gamma"]),
        colorbar_title = L"\mathbf{Q}",
        guidefontfamily="cm", fontfamily="cm",
        title = title)

    plot!(thetas*pi*2, repeat([4*pi/3], length(thetas)),
        c=cgrad(:darktest, rev=true)[data["triqs_data"]["path_ticks"][2]/data["triqs_data"]["path_plot"][end]], l=:dash, lw=2,
        label=L"K", linealpha=0.5)

    plot!(thetas*pi*2, repeat([2*pi/sqrt(3)], length(thetas)),
    c=cgrad(:darktest, rev=true)[data["triqs_data"]["path_ticks"][3]/data["triqs_data"]["path_plot"][end]], l=:dash, lw=2,
    label=L"M", linealpha=0.5)

    return p
end

########* Plotting peak strength *#########
function PlotJs(data::Dict, Js::Vector{Float64}, Qs::Vector{Vector{Float64}}, thetas::Vector{Float64} ;
    skip_every::Int64 = 1, title::AbstractString = "")

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
        markersize = 5, markerstrokealpha = 0.25,
        m=cgrad(:darktest, rev=true), zcolor = colors[1:skip_every:end], clims=(0, 1),
        colorbar_ticks = (data["triqs_data"]["path_ticks"]/data["triqs_data"]["path_plot"][end], [L"\Gamma", L"K", L"M", L"\Gamma"]),
        colorbar_title = L"\mathbf{Q}",
        guidefontfamily="cm", fontfamily="cm", title = title)

    return p
end

function PlotPhaseDiagram(t3::Float64, Bxs::Vector{Float64})


end

t3s = [-0.1, -0.2, -0.3, -0.33, -0.35, -0.37, -0.4, -0.5]

t3 = -0.2
Bx = -4.0


data = load("/home/anjishnubose/Research/Repos/RPA.jl/saves/data/bcao_Dirac_t3=$(round(t3, digits=2))_Bx=$(round(Bx, digits=2))_combined.jld2")
data = data["beta=20.0_mu=0.0"]

thetas = collect(LinRange(0, 1, 401))
labels = ["theta = $(round(theta, digits=3))*2*pi" for theta in thetas]

ks = data["triqs_data"]["contracted"]
ks = Vector{Float64}[eachrow(ks)...]

Js = Float64[]
Qs = Vector{Float64}[]

for label in labels
    push!(Js, data[label]["critical strength"])
    push!(Qs, data[label]["maximum momentum"])
end

pJ = PlotJs(data, Js, Qs, thetas)
title!(pJ, L"t_3 = %$(round(t3, digits=2)), B_x = %$(round(Bx, digits=2))")
savefig(pJ, "/home/anjishnubose/Research/Repos/RPA.jl/saves/plots/bcao_Dirac_t3=$(round(t3, digits=2))_Bx=$(round(Bx, digits=2))_Jcrit.pdf")
pQ = PlotQs(data, Qs, thetas)
title!(pQ, L"t_3 = %$(round(t3, digits=2)), B_x = %$(round(Bx, digits=2))")
savefig(pQ, "/home/anjishnubose/Research/Repos/RPA.jl/saves/plots/bcao_Dirac_t3=$(round(t3, digits=2))_Bx=$(round(Bx, digits=2))_Qcrit.pdf")


# for t3 in t3s
#     data = load("/home/anjishnubose/Research/Repos/RPA.jl/saves/data/bcao_Dirac_t3=$(round(t3, digits=2))_combined.jld2")
#     data = data["beta=20.0_mu=0.0"]

#     thetas = collect(LinRange(0, 1, 401))
#     # thetas = collect(LinRange(0.375, 0.425, 101))
#     labels = ["theta = $(round(theta, digits=3))*2*pi" for theta in thetas]

#     ks = data["triqs_data"]["contracted"]
#     ks = Vector{Float64}[eachrow(ks)...]

#     Js = Float64[]
#     Qs = Vector{Float64}[]

#     for label in labels
#         push!(Js, data[label]["critical strength"])
#         push!(Qs, data[label]["maximum momentum"])
#     end

#     pQ = PlotQs(data, Qs, thetas, title = L"t_3 = %$(round(t3, digits=2))")
#     pJ = PlotJs(data, Js, Qs, thetas, title = L"t_3 = %$(round(t3, digits=2))")

#     savefig(pQ, "/home/anjishnubose/Research/Repos/RPA.jl/saves/plots/bcao_Dirac_t3=$(round(t3, digits=2))XXZ_Qcrit.pdf")
#     savefig(pJ, "/home/anjishnubose/Research/Repos/RPA.jl/saves/plots/bcao_Dirac_t3=$(round(t3, digits=2))XXZ_Jcrit.pdf")

# end
