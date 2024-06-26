using ArgParse, JLD2, TightBindingToolkit, YAML

function parse_commandline()

    settings = ArgParseSettings()

    @add_arg_table settings begin
        "--input"
        help = "directory of the input yml file."
        arg_type = String
        required = true
    end

    return parse_args(settings)
end


if abspath(PROGRAM_FILE) == @__FILE__

    include("../RPA.jl")
    using .RPA

    parsed_args = parse_commandline()
    input = YAML.load_file(parsed_args["input"])

    model = load(input["unitcell"]["julia"])

    unitcell = model["unit cell"]
    parameters = model["parameters"]
    triqs_input = input["unitcell"]["julia"][1:end-5] * ".npz"

    parse_unitcell(unitcell, triqs_input)

    input["unitcell"]["triqs"] = triqs_input
    YAML.write_file(parsed_args["input"], input)

    command = `conda run -n $(input["triqs_environment"]) python ./Bare/run_bare.py $(parsed_args["input"])`
    run(command)

end
