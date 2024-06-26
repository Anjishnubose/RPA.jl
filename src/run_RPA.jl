using ArgParse, JLD2, TightBindingToolkit, YAML, NPZ

function parse_commandline()

    settings = ArgParseSettings()

    @add_arg_table settings begin
        "--input"
            help = "directory of the input yml file."
            arg_type = String
            required = true
        "--run_bare"
            help = "does the bare susceptibility calculation need to be run in TRIQS."
            arg_type = Bool
            default = false
    end

    return parse_args(settings)
end


if abspath(PROGRAM_FILE) == @__FILE__

    include("./RPA.jl")
    using .RPA

    parsed_args = parse_commandline()
    input = YAML.load_file(parsed_args["input"])

    localDim = length(input["directions"])

    model = load(input["unitcell"]["julia"])
    println("Unit Cell loaded!")

    unitcell = model["unit cell"]
    parameters = model["parameters"]
    subs = length(model["unit cell"].basis)

    interactions = load(input["interactions"])["parameters"]
    interaction_values = getproperty.(interactions, :value)
    n_interactions = length(interaction_values[begin])

    #####* running the bare
    if parsed_args["run_bare"]
        println("Running the bare susceptibility calculation in TRIQS...")

        command = `julia --project=../Project.toml ./Bare/run_bare.jl --input=$(parsed_args["input"])`
        run(command)

        println("Bare susceptibility calculation complete!")
    end

    for mu in input["mus"]["values"]
        println("Working on mu = $(mu)...")
        triqs_data = npzread(input["output"] * "_beta=$(input["beta"])_mu=$(round(mu, digits=3)).npz")

        primitives = dress_primitives(triqs_data)

        ks_contracted = triqs_data["contracted"]
        ks = Vector{eltype(ks_contracted)}[eachrow(ks_contracted)...]

        chis = combine_chis(triqs_data; directions = input["directions"], subs = subs)

        println("Starting RPA...")
        for ind in 1:n_interactions
            value = getindex.(interaction_values, ind)
            push!.(getproperty.(interactions, :value), value)
            lookup = Lookup(interactions)

            instability = find_instability(chis, ks; primitives = primitives,
                subs = subs, localDim=localDim,
                lookup = lookup)

            save(input["output"] * "_beta=$(input["beta"])_mu=$(round(mu, digits=3))_interactionID=$(ind).jld2",
                    instability)
        end

    end

end
