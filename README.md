# RPA
Codebase to calculate bare susceptiblity of a lattice model (using TRIQS), and then perform general random phase approximation (RPA) on it to get instabilities of the model.

## Installation
- Install the [triqs_tprf](https://triqs.github.io/tprf/2.1.x/install.html) package in a python venv (with default name "triqs"). Note: triqs_tprf does not work on Windows currently.
- Add the julia repo locally in Package mode.

## Usage
- Create the free Hamiltonian first. This is done using the [TightBindingToolkit.jl](https://github.com/Anjishnubose/TightBindingToolkit.jl) interface. Some examples are given in [./examples/models](https://github.com/Anjishnubose/RPA.jl/tree/main/examples/models).
- Define the interactions which are to be used for the RPA. This can also done in the same interface, with examples given in [./examples/interactions](https://github.com/Anjishnubose/RPA.jl/tree/main/examples/interactions).
- Write a parent input file with details of parameters required for the calculation like the inverse temperature, number of k-points, number of matsubara frequencies etc. Refer to [./Inputs](https://github.com/Anjishnubose/RPA.jl/tree/main/Inputs) for formatting.
- Run the following command
  ```julia
  julia --project=../Project.toml --heap-size-hint=4G --input="../Inputs/name_of_input.yml" --run_bare=true
  ```
