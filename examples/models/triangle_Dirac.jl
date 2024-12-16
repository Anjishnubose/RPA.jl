using Plots, LaTeXStrings, TightBindingToolkit, JLD2

"""
This script sets up a simple triangular lattice with the pi/0 flux state and a 3-sublattice spin order.
"""

"""
The lattice primitive vector for a triangular lattice : a1 and a2.
"""
const a1 = [ 1.0, 0.0 ]
const a2 = [ 1/2 , sqrt(3)/2 ]

const a1New = 2*a2
const a2New = a1
UC = UnitCell( [a1New , a2New] , 2) ##### localDim=2 since we are working with spin-1/2 particles now

"""
Unit cell has 2 sub-lattices.
"""
const b1 = [ 0.0 , 0.0 ]
const b2 = b1+a2
AddBasisSite!.( Ref(UC) , [b1, b2])

"""
Adding structure to the lattice now, through the bond objects.
"""
SpinVec     =   SpinMats(1//2) ##### Working with spin-1/2
################ Flux pi-0 hoppings #################
const t1  =   -1.0
const NNdistance  =   1.0
t1Param =  Param(t1, 2)

AddAnisotropicBond!(t1Param, UC, 1, 1, [ 0, 1],  -SpinVec[4], NNdistance, "dirac hoppings")
AddAnisotropicBond!(t1Param, UC, 1, 2, [ 0, 0],  +SpinVec[4], NNdistance, "dirac hoppings")
AddAnisotropicBond!(t1Param, UC, 1, 2, [ 0,-1],  -SpinVec[4], NNdistance, "dirac hoppings")

AddAnisotropicBond!(t1Param, UC, 2, 2, [ 0, 1],  +SpinVec[4], NNdistance, "dirac hoppings")
AddAnisotropicBond!(t1Param, UC, 2, 1, [ 1, 0],  +SpinVec[4], NNdistance, "dirac hoppings")
AddAnisotropicBond!(t1Param, UC, 2, 1, [ 1,-1],  +SpinVec[4], NNdistance, "dirac hoppings")
############### Zeeman field #################
const zeemanZ = 0.0
zeemanZParam = Param(zeemanZ, 2)
AddIsotropicBonds!(zeemanZParam, UC , 0.0, SpinVec[3] , "zeeman field along z")

params = [t1Param, zeemanZParam]
CreateUnitCell!(UC, [t1Param, zeemanZParam])

#####* Saving the unit cell in a JLD2 file
file_name = "../../saves/models/triangular_Dirac_Bz=$(round(zeemanZ, digits=2)).jld2"
save(file_name, Dict("unit cell" => UC, "parameters" => params))
