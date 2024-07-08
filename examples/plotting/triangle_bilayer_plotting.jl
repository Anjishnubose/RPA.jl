using JLD2, Plots, LaTeXStrings, LinearAlgebra

datas = load("/home/anjishnubose/Research/Repos/RPA.jl/saves/data/triangle_bilayer_combined.jld2")
const thetas = collect(LinRange(0, 1, 51))
labels = ["theta = $(round(theta/4, digits=3))*2*pi" for theta in thetas]

mus = Float64[]
fillings = Float64[]
keys = []

for (key, data) in datas

    push!(mus, data["mu"])
    push!(fillings, data["triqs_data"]["filling"])
    push!(keys, key)
end

args = sortperm(mus)
mus = mus[args]
fillings = fillings[args]
keys = keys[args]

const index = 7
mu = mus[index]
filling = fillings[index]
key = keys[index]

data = datas[key]
ks = data["triqs_data"]["contracted"]
ks = Vector{Float64}[eachrow(ks)...]


println("mu = $(mu), filling = $(filling)")

strengths = Float64[]
Qs = Vector{Float64}[]

for label in labels
    push!(strengths, data[label]["critical strength"])
    push!(Qs, data[label]["maximum momentum"])
end

const skip_every = 2

pyplot()
colors = data["triqs_data"]["path_plot"][indexin(Qs, ks)] / data["triqs_data"]["path_plot"][end]
plot(thetas[1:skip_every:end]*pi/2, norm.(Qs[1:skip_every:end]), label = L"|Q_{crit}|",
    framestyle=:box, gridalpha=0.6,
    xticks = ([0.0, 0.25*pi, 0.5*pi, 0.75*pi, pi, 1.25*pi, 1.5*pi, 1.75*pi], ["0", L"\frac{\pi}{4}", L"\frac{\pi}{2}", L"\frac{3\pi}{4}", L"\pi", L"\frac{5\pi}{4}", L"\frac{3\pi}{2}", L"\frac{7\pi}{4}"]),
    clabel = L"\theta",
    marker=:o, lw=2.0, proj=:polar, legend_position=:bottomleft,
    m=cgrad(:darktest, rev=true), zcolor = colors[1:skip_every:end], clims=(0, 1),
    colorbar_ticks = (data["triqs_data"]["path_ticks"]/data["triqs_data"]["path_plot"][end], [L"\Gamma", L"K", L"M", L"\Gamma"]),
    colorbar_title = L"Q",)

plot!(thetas*pi/2, repeat([4*pi/3], length(thetas)),
    c=cgrad(:darktest, rev=true)[data["triqs_data"]["path_ticks"][2]/data["triqs_data"]["path_plot"][end]], l=:dash, lw=2,
    label=L"K", linealpha=0.5)

plot!(thetas*pi/2, repeat([2*pi/sqrt(3)], length(thetas)),
c=cgrad(:darktest, rev=true)[data["triqs_data"]["path_ticks"][3]/data["triqs_data"]["path_plot"][end]], l=:dash, lw=2,
label=L"M", linealpha=0.5)

# plot(thetas[1:skip_every:end]*pi/2, strengths[1:skip_every:end], label = L"I_{crit}",
# framestyle=:box, gridalpha=0.6,
# xticks = ([0.0, 0.25*pi, 0.5*pi, 0.75*pi, pi, 1.25*pi, 1.5*pi, 1.75*pi], ["0", L"\frac{\pi}{4}", L"\frac{\pi}{2}", L"\frac{3\pi}{4}", L"\pi", L"\frac{5\pi}{4}", L"\frac{3\pi}{2}", L"\frac{7\pi}{4}"]),
# clabel = L"\theta",
# marker=:o, markersize=4.0, lw=2.0,
# proj=:polar,
# legend_position=:bottomleft,
# m=:darktest, zcolor = colors[1:skip_every:end], clims=(0, 1),
# colorbar_ticks = (data["triqs_data"]["path_ticks"]/data["triqs_data"]["path_plot"][end], [L"\Gamma", L"K", L"M", L"\Gamma"]),
# colorbar_title = L"Q",)

# plot!(thetas*pi/2, repeat([2*pi/sqrt(3)], length(thetas)),
#     c=:green, l=:dash, lw=2,
#     label=L"M", linealpha=0.75)
