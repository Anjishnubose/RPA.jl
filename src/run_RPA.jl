using ArgParse, JLD2, TightBindingToolkit, YAML, NPZ, LaTeXStrings, Plots

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
        "--plot_RPA"
            help = "does the RPA calculation need to be plotted."
            arg_type = Bool
            default = false
        "--save_individual"
            help = "does the individual RPA calculations need to be saved."
            arg_type = Bool
            default = false
    end

    return parse_args(settings)
end


if abspath(PROGRAM_FILE) == @__FILE__

    include("./RPA.jl")
    using .RPA

    parsed_args = parse_commandline()

    #####* reading the parent input file
    input = YAML.load_file(parsed_args["input"])
    localDim = length(input["directions"])
    model = load(input["unitcell"]["julia"])
    println("Unit Cell loaded!")

    #####* Extracting the unit cell and parameters from the model
    unitcell = model["unit cell"]
    parameters = model["parameters"]
    subs = length(model["unit cell"].basis)

    #####* Reading interaction data from file : parameters and different values to run RPA on
    interaction_data = load(input["interactions"])
    interaction_params = interaction_data["parameters"]
    interaction_values = interaction_data["values"]

    #####* running the bare susceptibility calculation in TRIQS if needed
    if parsed_args["run_bare"]
        println("Running the bare susceptibility calculation in TRIQS...")

        command = `julia --project=../Project.toml ./Bare/run_bare.jl --input=$(parsed_args["input"])`
        run(command)
        input = YAML.load_file(parsed_args["input"])
        command = `conda run -n $(input["triqs_environment"]) python ./Bare/plot_bare.py $(parsed_args["input"])`
        run(command)

        println("Bare susceptibility calculation complete!")
    end

    input = YAML.load_file(parsed_args["input"])
    #####* plot labels for k-space plots of bands and susceptibility
    path_labels = [L"%$(label)" for label in input["k_labels"]]
    push!(path_labels, path_labels[1])

    CombinedOutput = Dict()

    for mu in input["mus"]["values"]
        println("Working on mu = $(mu)...")
        parent_label = "beta=$(input["beta"])_mu=$(round(mu, digits=3))"
        triqs_data = npzread(input["output"] * "_$(parent_label).npz")

        CombinedOutput[parent_label] = Dict()
        CombinedOutput[parent_label]["mu"] = mu
        CombinedOutput[parent_label]["beta"] = input["beta"]
        CombinedOutput[parent_label]["triqs_data"] = triqs_data

        primitives = dress_primitives(triqs_data)

        ks_contracted = triqs_data["contracted"]
        ks = Vector{eltype(ks_contracted)}[eachrow(ks_contracted)...]
        #####* combining different chis to form the full susceptibility matrix
        chis = combine_chis(triqs_data; directions = input["directions"], subs = subs)
        CombinedOutput[parent_label]["combined_chis"] = chis

        println("Starting RPA...")
        for (label, value) in interaction_values
            #####* setting the interaction values
            push!.(getproperty.(interaction_params, :value), value)
            lookup = Lookup(interaction_params)

            CombinedOutput[parent_label]["$(label)_interaction"] = interaction(1.0, ks; primitives = primitives,
                                                                                subs = subs, localDim=localDim,
                                                                                lookup = lookup)

            instability = find_instability(chis, ks; primitives = primitives,
                subs = subs, localDim=localDim,
                lookup = lookup,
                upper=20.0,)

            critical = instability["critical strength"]
            strengths = collect(LinRange(0.0, 0.99*critical, 6))[2:end]

            if parsed_args["plot_RPA"]
                println("Plotting RPA for $(label)...")

                for strength in strengths
                    p = plot_chi(chis, strength, ks; primitives = primitives,
                        subs = subs, localDim=localDim,
                        lookup = lookup, path_plot = triqs_data["path_plot"],
                        path_ticks = triqs_data["path_ticks"],
                        path_labels = path_labels)

                    savefig(input["plots"] * "_chi_beta=$(input["beta"])_mu=$(round(mu, digits=3))_interactionID=$(label)_strength=$(round(strength, digits=3)).png")
                end
            end

            CombinedOutput[parent_label][label] = instability

            if parsed_args["save_individual"]
                save(input["output"] * "_beta=$(input["beta"])_mu=$(round(mu, digits=3))_interactionID=$(label).jld2",
                    instability)
            end
        end
        println("RPA complete!")

    end

    save(input["output"] * "_combined.jld2", CombinedOutput)

end
