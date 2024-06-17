module parse_model

using TightBindingToolkit, LinearAlgebra, NPZ

export parse_unitcell


function parse_unitcell(uc::UnitCell ;)::Dict{String, Any}

    output = Dict{String, Any}()

    primitives = [[uc.primitives[1] ; 0.0], [uc.primitives[2] ; 0.0]]
    output["units"] = reduce(hcat, primitives)

    positions = repeat(uc.basis, inner = uc.localDim)
    output["orbital_positions"] = reduce(hcat, positions)

    lookup = Lookup(uc)
    n  = uc.localDim * length(uc.basis)
    hoppings = Dict()

    for (key, value) in lookup
        i, j, offset = key
        offset = Tuple(offset)

        b1  =   uc.localDim * (i - 1) + 1
        b2  =   uc.localDim * (j - 1) + 1

        if !haskey(hoppings, offset)
            hoppings[offset] = zeros(ComplexF64, n, n)
            hoppings[.-(offset)]= zeros(ComplexF64, n, n)
        end

        hoppings[offset][b1 : b1 + uc.localDim - 1, b2 : b2 + uc.localDim - 1] += value
        hoppings[.-(offset)][b2 : b2 + uc.localDim - 1, b1 : b1 + uc.localDim - 1] += collect(value')
    end

    offsets = collect(keys(hoppings))
    hopping_matrices = zeros(ComplexF64, length(offsets), n, n)
    for (i, offset) in enumerate(offsets)
        hopping_matrices[i, :, :] = hoppings[offset]
    end

    output["hopping matrices"] = hopping_matrices
    output["hopping offsets"] = reduce(hcat, collect.(offsets))

    return output
end


function parse_unitcell(uc::UnitCell , fileName::String)

    output = parse_unitcell(uc)
    npzwrite(fileName, output)

end


end
