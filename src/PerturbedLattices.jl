module PerturbedLattices

using LinearAlgebra
using Random
using Distributions
using StaticArrays
using Plots
using Optim
using JLD2
using Infinity

# Core types
export PerturbedLatticeModel

# Main functions
export create_grid, iterate!, shift!
export local_energy, local_energy_vec, local_energy_MC
export estimate_int_DLR, DLR, plot_DLR
export minimize_DLR_squared, minimize_DLR_squared_nelder
export points_in_window, create_grid_for_est_int
export DLR_W_n, minimize_DLR_squared_nelder_W_n
export plot_points, plot_point_grid_connection

# Include submodules
include("types.jl")
include("grid.jl")
include("energy.jl")
include("simulation.jl")
include("windowed_estimation.jl")
include("visualization.jl")

end 