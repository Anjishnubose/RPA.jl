using JLD2, Plots, LaTeXStrings, LinearAlgebra

datas = load("/home/anjishnubose/Research/Repos/RPA.jl/saves/data/triangle_NN_combined.jld2")
label = "NN_repulsive_density-density"

mus = Float64[]
fillings = Float64[]

Vs = Float64[]
Qs = Vector{Float64}[]

for (key, data) in datas

    push!(mus, data["mu"])
    push!(fillings, data["triqs_data"]["filling"])

    ks = data["triqs_data"]["contracted"]
    ks = Vector{Float64}[eachrow(ks)...]

    push!(Vs, data[label]["critical strength"])
    push!(Qs, data[label]["maximum momentum"])

end

args = sortperm(mus)
mus = mus[args]
fillings = fillings[args]
Vs = Vs[args]
Qs = Qs[args]
