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
jh = get!(params, "jh", 2.0)
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
##Functions that will be useful for adding anisotropic bonds
weiss_neel(v) = [sin(pi * (norm(v) / (SkXSize))) * v[1] / norm(v), sin(pi * (norm(v) / (SkXSize))) * v[2] / norm(v), cos(pi * (norm(v) / (SkXSize)))]
weiss_bloch(v) = [sin(pi * (norm(v) / (SkXSize))) * v[2] / norm(v), sin(pi * (norm(v) / (SkXSize))) * -v[1] / norm(v), cos(pi * (norm(v) / (SkXSize)))]
weiss = Dict("Neel" => weiss_neel, "Bloch" => weiss_bloch)
sigmav(i, j) = 2 .* [su2spin[1][i, j], su2spin[2][i, j], su2spin[3][i, j]]
s11 = sigmav(1, 1)
s12 = sigmav(1, 2)
s21 = sigmav(2, 1)
s22 = sigmav(2, 2)
intermat(s) = [dot(s, s11) dot(s, s12); dot(s, s21) dot(s, s22)]

for (ind, bas) in enumerate(UC.basis)
    closest = [bas, bas - a1, bas - a2, bas - a1 - a2, bas + a1, bas + a2, bas + a1 + a2, bas + a1 - a2, bas - a1 + a2]
    minimal = findmin(x -> norm(x), closest)[2]
    if (SkXSize - 1) < norm(closest[minimal]) < SkXSize
        mat = intermat(normalize(weiss[SkX](closest[minimal]) + weiss[SkX](-closest[minimal])))
    else
        spn = weiss[SkX](closest[minimal])
        replace!(spn, NaN => 0.0)
        mat = intermat(normalize(spn))
    end
    AddAnisotropicBond!(jhParam, UC, ind, ind, [0, 0], mat, 0.0, "Hunds")
end

#####* Adding the parameter to the unit cell
CreateUnitCell!(UC, HoppingParams)


#####* Saving the unit cell in a JLD2 file
# file_name = "/scratch/a/aparamek/andykh/Data/Monolayer_Data/RPA/SkX_4.jld2"
# change to scratch
#save(file_name, Dict("unit cell" => UC, "parameters" => HoppingParams))
kSize = 10*6 + 3
bz = BZ(kSize)
FillBZ!(bz, UC)
H = Hamiltonian(UC, bz)
DiagonalizeHamiltonian!(H)
Mdl = Model(UC, bz, H; filling=filling)
SolveModel!(Mdl; get_gap=true)

bands = Plot_Band_Structure!(Mdl, [bz.HighSymPoints["G"], bz.HighSymPoints["K1"], bz.HighSymPoints["M2"]], labels=["G", "K1", "M2"], plot_legend=false)
