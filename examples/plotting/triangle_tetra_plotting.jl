using JLD2, Plots, LaTeXStrings, LinearAlgebra

include("/home/anjishnubose/Research/Repos/RPA.jl/src/RPA.jl")
using .RPA


datas = load("/home/anjishnubose/Research/Repos/RPA.jl/saves/data/triangle_tetra_combined.jld2")
label = "NN_repulsive_density-density"

mus = Float64[]
fillings = Float64[]

Vs = Float64[]
Qs = Vector{Float64}[]

key1 = collect(keys(datas))[begin]
data1 = datas[key1]
primitives = dress_primitives(data1["triqs_data"])
primitives = [[primitives[1]; 0.0], [primitives[2]; 0.0], [0.0, 0.0, 1.0]]
reciprocal = dress_reciprocal(data1["triqs_data"])

ks = data1["triqs_data"]["contracted"]
ks = Vector{Float64}[eachrow(ks)...]

for (key, data) in datas

    push!(mus, data["mu"])
    push!(fillings, data["triqs_data"]["filling"])

    push!(Vs, data[label]["critical strength"])
    push!(Qs, data[label]["maximum momentum"])

end

args = sortperm(mus)
mus = mus[args]
fillings = fillings[args]
Vs = Vs[args]
Qs = Qs[args]
