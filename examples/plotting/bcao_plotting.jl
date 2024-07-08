using JLD2, Plots, LaTeXStrings, LinearAlgebra

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

# plot(thetas .* (2*pi), Js, proj=:polar, gridalpha = 0.75, lw=2.0, label=L"J_{crit}", marker=:circle, title=L"t_1=%$(-1.0)\,,t_3=%$(-0.2), \beta=%$(20.0)", m=:darktest, zcolor = indexin(Qs, ks) ./ 300, clims=(0, 1))


pyplot()

const skip_every = 1
colors = data["triqs_data"]["path_plot"][indexin(Qs, ks)] / data["triqs_data"]["path_plot"][end]

# plot(thetas[1:skip_every:end]*pi*2, norm.(Qs[1:skip_every:end]), label = L"|Q_{crit}|",
#     framestyle=:box, gridalpha=0.6,
#     background_color_legend=RGBA(1.0, 1.0, 1.0, 0.9),
#     xticks = ([0.0, 0.25*pi, 0.5*pi, 0.75*pi, pi, 1.25*pi, 1.5*pi, 1.75*pi], ["0", L"\frac{\pi}{4}", L"\frac{\pi}{2}", L"\frac{3\pi}{4}", L"\pi", L"\frac{5\pi}{4}", L"\frac{3\pi}{2}", L"\frac{7\pi}{4}"]),
#     clabel = L"\theta",
#     marker=:o, lw=2.0, proj=:polar, legend_position=:bottomleft,
#     m=cgrad(:darktest, rev=true), zcolor = colors[1:skip_every:end], clims=(0, 1),
#     colorbar_ticks = (data["triqs_data"]["path_ticks"]/data["triqs_data"]["path_plot"][end], [L"\Gamma", L"K", L"M", L"\Gamma"]),
#     colorbar_title = L"\mathbf{Q}",)

# plot!(thetas*pi*2, repeat([4*pi/3], length(thetas)),
#     c=cgrad(:darktest, rev=true)[data["triqs_data"]["path_ticks"][2]/data["triqs_data"]["path_plot"][end]], l=:dash, lw=2,
#     label=L"K", linealpha=0.5)

# plot!(thetas*pi*2, repeat([2*pi/sqrt(3)], length(thetas)),
# c=cgrad(:darktest, rev=true)[data["triqs_data"]["path_ticks"][3]/data["triqs_data"]["path_plot"][end]], l=:dash, lw=2,
# label=L"M", linealpha=0.5)

plot(thetas[1:skip_every:end]*pi*2, (Js[1:skip_every:end]), label = L"|J_{crit}|",
    framestyle=:box, gridalpha=0.6,
    background_color_legend = RGBA(1.0, 1.0, 1.0, 0.4),
    xticks = ([0.0, 0.25*pi, 0.5*pi, 0.75*pi, pi, 1.25*pi, 1.5*pi, 1.75*pi], ["0", L"\frac{\pi}{4}", L"\frac{\pi}{2}", L"\frac{3\pi}{4}", L"\pi", L"\frac{5\pi}{4}", L"\frac{3\pi}{2}", L"\frac{7\pi}{4}"]),
    clabel = L"\theta",
    marker=:o, lw=2.0, proj=:polar, legend_position=:topright,
    m=cgrad(:darktest, rev=true), zcolor = colors[1:skip_every:end], clims=(0, 1),
    colorbar_ticks = (data["triqs_data"]["path_ticks"]/data["triqs_data"]["path_plot"][end], [L"\Gamma", L"K", L"M", L"\Gamma"]),
    colorbar_title = L"\mathbf{Q}",)
