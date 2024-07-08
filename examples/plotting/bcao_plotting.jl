using JLD2, Plots, LaTeXStrings

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

plot(thetas .* (2*pi), Js, proj=:polar, gridalpha = 0.75, lw=2.0, label=L"J_{crit}", marker=:circle, title=L"t_1=%$(-1.0)\,,t_3=%$(-0.2), \beta=%$(20.0)", m=:darktest, zcolor = indexin(Qs, ks) ./ 300, clims=(0, 1))
