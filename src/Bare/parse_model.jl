module parse_model

using JLD2, TightBindingToolkit, LinearAlgebra, YAML

export parse_unitcell

function parse_unitcell(uc::UnitCell ;
        orbital_labels::Vector{String} = Ref("orb = ") .* string.(1:uc.localDim))::Dict

    output = Dict{String, Any}()
    output["units"] = [[uc.primitives[1] ; 0.0], [uc.primitives[2] ; 0.0]]
    output["orbital_positions"] = repeat(uc.basis, inner = uc.localDim)

    sublattice_names = Ref("sub = ") .* string.(1:length(uc.basis))
    orbital_names = [sublattice * "_" * orbital for orbital in orbital_labels, sublattice in sublattice_names]

    output["orbital_names"] = reshape(orbital_names, length(orbital_names))

    lookup = Lookup(uc)
    n  = uc.localDim * length(uc.basis)
    hoppings = Dict{Tuple{Int, Int}, Matrix{ComplexF64}}()

    for (key, value) in lookup
        i, j, offset = key
        offset = Tuple(offset)

        b1  =   uc.localDim * (i - 1) + 1
        b2  =   uc.localDim * (j - 1) + 1

        if !haskey(hoppings, offset)
            hoppings[offset] = zeros(ComplexF64, n, n)
            hoppings[.-(offset)] = zeros(ComplexF64, n, n)
        end

        hoppings[offset][b1 : b1 + uc.localDim - 1, b2 : b2 + uc.localDim - 1] += value
        hoppings[.-(offset)][b2 : b2 + uc.localDim - 1, b1 : b1 + uc.localDim - 1] += collect(value')

    end

    output["hoppings"] = hoppings

    return output
end

function parse_unitcell(uc::UnitCell , fileName::String;
    orbital_labels::Vector{String} = Ref("orb = ") .* string.(1:uc.localDim))::Dict

    output = parse_unitcell(uc, orbital_labels = orbital_labels)
    YAML.write_file(fileName, output)

end





end
