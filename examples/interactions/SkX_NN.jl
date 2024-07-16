using JLD2, TightBindingToolkit, LinearAlgebra

#####* Define the honeycomb lattice unit vectors
params = Dict()
SkXSize = get!(params, "SkXSize", 2)
SkX = get!(params, "SkX", "Bloch")
SkX = "Bloch"
a1 = SkXSize / 2 * [-3.0, sqrt(3)]
a2 = SkXSize / 2 * [3.0, sqrt(3)]
l1 = [1.0, 0]
l2 = [-0.5, sqrt(3) / 2]
UC = UnitCell([a1, a2], 2, 2)
##Parameters
n = get!(params, "n", 20)
kSize = 6 * n + 3
t = get!(params, "t", 1.0)
jh = get!(params, "jh", 5.0)
U = get!(params, "U", 0.0)
##### Thermodynamic parameters
filling = get!(params, "filling", 12.5/24)
T = get!(params, "T", 0.0)
t1 = -t
t1Param = Param(t1, 2)
jhParam = Param(jh, 2)
HoppingParams = [t1Param, jhParam]
su2spin = SpinMats(1 // 2)

##Adding inner-hexagon structure
for j = 0:(SkXSize-1)
    for i = 0:(SkXSize*3-1)
        AddBasisSite!(UC, i .* l1 + j .* l2)
    end
end
AddIsotropicBonds!(t1Param, UC, 1.0, su2spin[4], "t1", checkOffsetRange=1)


function get_distances(primitives::Vector{Vector{Float64}}, check_range::Int64)

    distances = [norm(sum([i, j] .* primitives)) for i in -check_range:check_range, j in -check_range:check_range]
    distances = reshape(distances, length(distances))
    sort!(distances)

    bonds = Int64[1]
    unique_distances = [distances[1]]
    k=1
    for (d, distance) in enumerate(distances[2:end])
        if abs(distance-distances[d])<1e-6
            bonds[k] +=1
        else
            push!(bonds, 1)
            push!(unique_distances, distance)
            k+=1
        end
    end

    return (unique_distances, bonds)
end

distances, counts = get_distances([a1, a2], 10)

#####* Adding a parameter which tracks the nearest neighbor Heisenberg interaction in the spin-spin basis

VParams = Param{2, Float64}[]
for (d, distance) in enumerate(distances[2:end])

    if counts[d+1] % 6 ==0
        param   =   Param(4.0, 2)
        AddIsotropicBonds!(param, UC , distance, [1.0;;] / distance , "V_$(d)" ; checkOffsetRange = 10)
        push!(VParams, param)
    end

end


values = Dict()

params = VParams
values["Coulomb_Repulsion"] = repeat([1.0], length(params))

#####* Saving the unit cell in a JLD2 file
# file_name = "/scratch/a/aparamek/andykh/Data/Monolayer_Data/RPA/SkX_NN.jld2"
#save(file_name, Dict("parameters" => params))
# save(file_name, Dict("parameters" => params, "values" => values))
