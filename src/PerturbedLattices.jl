module PerturbedLattices

using LinearAlgebra
using Random
using Distributions
using StaticArrays
using Plots
using Optim
using JLD2
using Infinity

export PerturbedLatticeModel
export Grid
export StraussHamiltonian, HardCoreHamiltonian
export GaussianMoveModel, UniformMoveModel
export rand
# Core types
export PerturbedLatticeV1

# Main functions
export create_grid, iterate!, shift!
export local_energy
export points_in_window 
export DLR_W, fit
export plot_points, plot_point_grid_connection

const Point= Vector{Float64}
abstract type AbstractLattice end
abstract type AbstractPerturbedLatticeModel end

# Include submodules
include("grid.jl")
include("hamiltonian.jl")
include("move_model.jl")
#include("simulation.jl")
#include("visualization.jl")
#include("estimation.jl")
#include("data_creation.jl")
include("perturbed_lattice_model.jl")


include("v1/perturbedLatticeModel.jl")
include("v1/grid.jl")
include("v1/energy.jl")
include("v1/simulation.jl")
include("v1/visualization.jl")
include("v1/estimation.jl")
include("v1/data_creation.jl")

end