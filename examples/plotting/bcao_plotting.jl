using JLD2, Plots, LaTeXStrings, LinearAlgebra
include("../../src/RPA.jl")
using .RPA
gr()
#####################################
const cg = cgrad(:darktest, rev=true)
const color_ZZ = cg[0.64]
const color_FM = cg[0.0]
const color_SpiralK = cg[0.2]
const color_SpiralM = cg[0.8]
#####################################

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

function PlotPhaseDiagram(t3::Float64, Bxs::Vector{Float64}, thetas::Vector{Float64}, J0::Float64, skip_every::Int64 = 1)

    labels = ["theta = $(round(theta, digits=3))*2*pi" for theta in thetas]
    p = plot(framestyle=:box, grid=false,
                guidefont = "Computer Modern", tickfont = "Computer Modern", legendfont = "Computer Modern", titlefont = "Computer Modern",
                guidefontsize = 14, tickfontsize = 12, legendfontsize = 12, titlefontsize = 14,
                xlabel = L"J_3/J_1", ylabel = L"B_x/t_1", title = L"t_3 = %$(round(t3, digits=2))\,,J = %$(round(J0, digits=2))",
                )

    for Bx in Bxs
        fileName = "/home/anjishnubose/Research/Repos/RPA.jl/saves/data/bcao_Dirac_t3=$(round(t3, digits=2))_Bx=$(round(Bx, digits=2))_combined.jld2"
        data = load(fileName)
        data = data["beta=20.0_mu=0.0"]

        ks = data["triqs_data"]["contracted"]
        ks = Vector{Float64}[eachrow(ks)...]

        Js = Float64[]
        Qs = Vector{Float64}[]

        for label in labels
            push!(Js, data[label]["critical strength"])
            push!(Qs, data[label]["maximum momentum"])
        end
        m = cgrad(:darktest, rev=true)
        colors = []
        zcolors = data["triqs_data"]["path_plot"][indexin(Qs, ks)] / data["triqs_data"]["path_plot"][end]
        for (j, J) in enumerate(Js)
            if (J>J0)
                push!(colors, :white)
            else
                push!(colors, m[zcolors[j]])
            end
        end

        scatter!(tan.(thetas[1:skip_every:end]*2*pi), repeat([abs(Bx)], length(thetas))[1:skip_every:end],
                markersize=8, markerstrokealpha=0.1, markercoloralpha = 0.75, markerstrokewidth=0.5,
                markercolor = colors[1:skip_every:end], label = "")
    end


    return p
end

function relPhase(v::Vector{ComplexF64})
    if abs(v[1]) > 1e-2
        return angle(v[4]/v[1])
    else
        return angle(v[5]/v[2])
    end
end

function readData(Parentdata::Dict, thetas::Vector{Float64})
    data = Parentdata["beta=20.0_mu=0.0"]
    labels = ["theta = $(round(theta, digits=3))*2*pi" for theta in thetas]

    ks = data["triqs_data"]["contracted"]
    ks = Vector{Float64}[eachrow(ks)...]

    Js = Float64[]
    Qs = Vector{Float64}[]
    eigenvectors = Vector{ComplexF64}[]

    for label in labels
        push!(Js, data[label]["critical strength"])
        push!(Qs, data[label]["maximum momentum"])
        push!(eigenvectors, data[label]["maximum eigenvector"])
    end

    ks = data["triqs_data"]["contracted"]
    ks = Vector{Float64}[eachrow(ks)...]
    colors = data["triqs_data"]["path_plot"][indexin(Qs, ks)] / data["triqs_data"]["path_plot"][end]

    return Dict("Js" => Js, "Qs" => Qs, "eigenvectors" => eigenvectors, "colors" => colors, "ks" => ks)
end

function annotatewithbox!(
    fig::Plots.Plot,
    text::Plots.PlotText,
    x::Real, y::Real, Δx::Real, Δy::Real = Δx;
    color = :white,
    kwargs...)

    box = Plots.Shape(:rect)
    Plots.scale!(box, 1/sqrt(2), 1/sqrt(2))
    Plots.scale!(box, Δx, Δy)
    Plots.translate!(box, x, y)

    Plots.plot!(fig, box, c = color, linestroke = :black, label = false; kwargs...)
    Plots.annotate!(fig, x, y, text)

end


function parabola(x0::Float64, y0::Float64, x1::Float64)

    function p(x::Float64)
        return y0*(1-((x-x0)^2/(x1-x0)^2))
    end
    return p
end


function findAlpha(thetas::Vector{Float64}, Js::Vector{Float64}, chi::Float64=0.13, ar::Float64=1.0)
    alphas = (1 .+ chi*(Js ./ abs.(sec.(thetas)))) .^ (-1)
    J3s = ((1 .- alphas) ./ (1 .- ar*alphas)) .* tan.(thetas)
    return J3s, alphas
end


function PlotAlphaPD(data::Dict, thetas::Vector{Float64};
        starting::Int64 = 91, skip_every::Int64 = 1, ending::Int64 = length(thetas),
        chi::Float64 = 0.13,
        ar::Float64 = 1.0,
        classical_boundaries::Vector{Float64} = [-0.4, -0.25])

    output = readData(data, thetas)
    data = data["beta=20.0_mu=0.0"]
    Js = output["Js"][starting:skip_every:ending]
    Qs = output["Qs"][starting:skip_every:ending]

    xs, ys = findAlpha(thetas[starting:skip_every:ending]*2*pi, Js, chi, ar)

    spiral_start = findfirst(x -> abs(x-0.63) > 0.05 , output["colors"][starting:skip_every:ending])
    spiral_end = findfirst(x -> abs(x) < 0.05 , output["colors"][starting:skip_every:ending])
    spiral_start = (xs[spiral_start], ys[spiral_start])
    println(spiral_start)
    spiral_end = (xs[spiral_end-1], ys[spiral_end-1])

    shift = 0.0
    boundary1xs = collect(LinRange(spiral_start[1], classical_boundaries[1]+shift, 101))
    boundary1ys = parabola(spiral_start[1], spiral_start[2], classical_boundaries[1]+shift).(boundary1xs)

    shift = 0.01
    boundary2xs = collect(LinRange(spiral_end[1], classical_boundaries[2]+shift, 101))
    boundary2ys = parabola(spiral_end[1], spiral_end[2], classical_boundaries[2]+shift).(boundary2xs)
    gr()
    p = scatter(
    ylabel = L"\alpha",
    framestyle=:box, grid = false, label = "",
    xlabel = L"J_3/J_1",
    guidefont = "Computer Modern", legendfont = "Computer Modern", tickfont = "Computer Modern",
    guidefontsize = 14, tickfontsize = 12, legendfontsize = 12,
    ylims = (-0.15, 1.05),
    xlims = (-1.05, 0.025))

    plot!(boundary1xs, boundary1ys, lw=4, c=color_ZZ, label="", linealpha = 0.6)
    plot!(boundary2xs, boundary2ys, lw=4, c=color_FM, label="", linealpha = 0.6)

    scatter!(xs, ys, label = L"\alpha_3/\alpha_1=%$(round(ar, digits=2))",
    marker=:circle, lw=2.0, legend_position=:topleft,
    markersize = 6, markerstrokealpha = 0.25,
    markercolors = cgrad(:darktest, rev=true)[output["colors"][starting:skip_every:ending]],
    ylims = (-0.15, 1.05),
    xlims = (-1.05, 0.025))

    zz_center = ((-1.05+classical_boundaries[1])/2, 0.0)
    zz_dimensions = (classical_boundaries[1]+1.05, 0.1)
    annotatewithbox!(p, text(L"Zig-Zag", :black, :center, 12, "Computer Modern"),
                zz_center..., zz_dimensions..., color = color_ZZ)

    ferro_center = ((0.025+classical_boundaries[2])/2, 0.0)
    ferro_dimensions = (classical_boundaries[2]-0.025, 0.1)
    annotatewithbox!(p, text(L"FM", :black, :center, 12, "Computer Modern"),
                ferro_center..., ferro_dimensions..., color = color_FM)

    spiral_center = ((classical_boundaries[1]+classical_boundaries[2])/2, 0.0)
    spiral_dimensions = (classical_boundaries[2]-classical_boundaries[1], 0.1)
    annotatewithbox!(p, text(L"Spiral", :black, :center, 12, "Computer Modern"),
                spiral_center..., spiral_dimensions..., color = color_SpiralM)

    annotate!(-0.5, -0.1, text("Classical", :black, :center, 12, "Computer Modern"))
    annotate!(-0.5, 0.95, text("Dirac SL", :black, :center, 12, "Computer Modern"))

    return p
end



const t3 = -0.2
const Bx = 0.0

# data = load("./saves/data/bcao_Dirac_t3=$(round(t3, digits=2))_Bx=$(round(Bx, digits=2))_combined.jld2")
data = load("./saves/data/bcao_Dirac_t3=$(round(t3, digits=2))_combined.jld2")
thetas = collect(LinRange(0.25, 0.5, 101))
labels = ["theta = $(round(theta, digits=3))*2*pi" for theta in thetas]

# data = readData(data, thetas)

# const skip_every = 1
# title = ""
# pyplot()
# p = plot(thetas[1:skip_every:end], abs.(relPhase.(eigenvectors[1:skip_every:end])/pi), ylabel = L"\phi/\pi",
#     framestyle=:box, grid = false, label = "",
#     background_color_legend = RGBA(1.0, 1.0, 1.0, 0.4),
#     xlabel = L"\theta/2\pi",
#     markershape=:circle, lw=2.0, legend_position=:bottomleft,
#     markersize = 6, markerstrokealpha = 0.25,
#     m=cgrad(:darktest, rev=true), zcolor = colors[1:skip_every:end], clims=(0, 1),
#     colorbar_ticks = (data["triqs_data"]["path_ticks"]/data["triqs_data"]["path_plot"][end], [L"\Gamma", L"K", L"M", L"\Gamma"]),
#     colorbar_title = L"\mathbf{Q}",
#     colorbar_titlefontsize=14,
#     colorbar_tickfontsize = 14,
#     # guidefont="Computer Modern", tickfont="Computer Modern", legendfont = "Computer Modern",
#     guidefontfamily = "cm", fontfamily = "cm", tickfontfamily = "cm", legendfontfamily = "cm",
#     guidefontsize = 14, tickfontsize = 12, legendfontsize = 12, titlefontsize = 14,
#     title = title,
#     ylims = (-0.05, 1.05))

#     savefig(p, "./bcao_Dirac_t3=$(round(t3, digits=2))_Bx=$(round(Bx, digits=2))_Pcrit.pdf")





# x1s = collect(LinRange(-0.4, -0.7, 101))
# x2s = collect(LinRange(-0.25, -0.65, 101))
# y1s = -(0.68/0.3)*(x1s .+ 0.4)
# y2s = -(0.68/0.4)*(x2s .+ 0.25)

# x1s = collect(LinRange(-0.4, -0.56, 101))
# x2s = collect(LinRange(-0.325, -0.52, 101))
# y1s = -(0.69/0.16)*(x1s .+ 0.4)
# y2s = -(0.69/0.195)*(x2s .+ 0.325)

# plot!(x1s, y1s, ls=:dash, lw=2, c=cgrad(:darktest, rev=true)[colors[101]], label="")
# plot!(x2s, y2s, ls=:dash, lw=2, c=cgrad(:darktest, rev=true)[colors[end]], label="")

# plot!(tan.(thetas[101:152]*2*pi), repeat([0.0], 52), markercolor = cgrad(:darktest, rev=true)[colors[101]],
#         markershape=:square, label="", markersize = 6, markerstrokealpha = 0.25,)
# plot!(tan.(thetas[153:170]*2*pi), repeat([0.0], 18), markercolor = cgrad(:darktest, rev=true)[colors[126]],
# markershape=:square, label="", markersize = 6, markerstrokealpha = 0.25,)
# plot!(tan.(thetas[171:201]*2*pi), repeat([0.0], 31), markercolor = cgrad(:darktest, rev=true)[colors[end]],
#         markershape=:square, label="", markersize = 6, markerstrokealpha = 0.25,)

# plot!(tan.(thetas[101:152]*2*pi), repeat([0.0], 52), markercolor = cgrad(:darktest, rev=true)[colors[101]],
# markershape=:square, label="", markersize = 6, markerstrokealpha = 0.25,)
# plot!(tan.(thetas[153:161]*2*pi), repeat([0.0], 9), markercolor = cgrad(:darktest, rev=true)[colors[139]],
# markershape=:square, label="", markersize = 6, markerstrokealpha = 0.25,)
# plot!(tan.(thetas[162:201]*2*pi), repeat([0.0], 40), markercolor = cgrad(:darktest, rev=true)[colors[end]],
# markershape=:square, label="", markersize = 6, markerstrokealpha = 0.25,)

# annotate!(-0.9, 0.05, text("Zig-Zag", :black, :center, 12))
# annotate!(-0.1, 0.05, text("FM", :black, :center, 12))
# annotate!(-0.375, -0.05, text("Spiral", :black, :center, 12))
# annotate!(-0.7, 0.85, text("Dirac SL", :black, :center, 12))



# savefig(p, "./bcao_Dirac_t3=$(round(t3, digits=2))_Bx=$(round(Bx, digits=2))_Jcrit.pdf")
# Bx = -1.0
# data = load("./saves/data/bcao_Dirac_t3=$(round(t3, digits=2))_Bx=$(round(Bx, digits=2))_combined.jld2")
# data = data["beta=20.0_mu=0.0"]

# ks = data["triqs_data"]["contracted"]
# ks = Vector{Float64}[eachrow(ks)...]

# Js = Float64[]
# Qs = Vector{Float64}[]

# for label in labels
#     push!(Js, data[label]["critical strength"])
#     push!(Qs, data[label]["maximum momentum"])
# end

# ks = data["triqs_data"]["contracted"]
# ks = Vector{Float64}[eachrow(ks)...]

# colors = data["triqs_data"]["path_plot"][indexin(Qs, ks)] / data["triqs_data"]["path_plot"][end]

# plot!(thetas[1:skip_every:end], (Js[1:skip_every:end]), ylabel = L"|J_{crit}|",
#     framestyle=:box, grid = false, label = L"B_x=1",
#     background_color_legend = RGBA(1.0, 1.0, 1.0, 0.4),
#     xlabel = L"\theta/2\pi",
#     markershape=:square, lw=2.0, legend_position=:bottomleft,
#     markersize = 6, markerstrokealpha = 0.25,
#     m=cgrad(:darktest, rev=true), zcolor = colors[1:skip_every:end], clims=(0, 1),
#     colorbar_ticks = (data["triqs_data"]["path_ticks"]/data["triqs_data"]["path_plot"][end], [L"\Gamma", L"K", L"M", L"\Gamma"]),
#     colorbar_title = L"\mathbf{Q}",
#     colorbar_tickfontsize = 14,
#     # guidefont="Computer Modern", tickfont="Computer Modern", legendfont = "Computer Modern",
#     guidefontfamily = "cm", fontfamily = "cm", tickfontfamily = "cm", legendfontfamily = "cm",
#     guidefontsize = 14, tickfontsize = 12, legendfontsize = 12, titlefontsize = 14,
#     title = title)

# pJ = PlotJs(data, Js, Qs, thetas)
# title!(pJ, L"t_3 = %$(round(t3, digits=2)), B_x = %$(round(Bx, digits=2))")
# savefig(pJ, "/home/anjishnubose/Research/Repos/RPA.jl/saves/plots/bcao_Dirac_t3=$(round(t3, digits=2))_Bx=$(round(Bx, digits=2))_Jcrit.pdf")
# pQ = PlotQs(data, Qs, thetas)
# title!(pQ, L"t_3 = %$(round(t3, digits=2)), B_x = %$(round(Bx, digits=2))")
# savefig(pQ, "/home/anjishnubose/Research/Repos/RPA.jl/saves/plots/bcao_Dirac_t3=$(round(t3, digits=2))_Bx=$(round(Bx, digits=2))_Qcrit.pdf")


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
