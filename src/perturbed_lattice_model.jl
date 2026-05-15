"""
    PerturbedLattice

A structure for managing perturbed lattice point processes in 2D or 3D.

# Fields
- `N::Int`: Half-size of the lattice (total points: (2N+1)^d)
- `d::Int`: Dimension (2 or 3)
- `NMC::Int`: Number of Monte Carlo iterations
- `H::Union{String, Nothing}`: Hamiltonian type ("S", "SHC", "HC", or nothing)
- `RS::Union{Float64, Nothing}`: Radius of Strauss potential
- `RHC::Union{Float64, Nothing}`: Radius of hardcore interaction
- `beta::Union{Float64, Nothing}`: Inverse temperature parameter
- `Distrib::Union{String, Nothing}`: Perturbation distribution ("Gauss", "Unif", or nothing)
- `BorneUnif::Union{Matrix{Float64}, Nothing}`: Bounds for uniform distribution
- `Cov::Union{Matrix{Float64}, Nothing}`: Covariance matrix for Gaussian distribution
- `grid::Vector{Vector{Float64}}`: Regular lattice points
- `points::Vector{Vector{Float64}}`: Perturbed lattice points
- `adjacencyS::Matrix{Int}`: Adjacency matrix for Strauss potential
- `adjacencyHC::Matrix{Int}`: Adjacency matrix for hardcore interaction
- `seed::Union{Int, Nothing}`: Random seed for reproducibility
- `rng::AbstractRNG`: Random number generator

"""

mutable struct PerturbedLatticeModel <: AbstractPerturbedLatticeModel
    h::AbstractHamiltonian

    move::AbstractMoveModel

    grid::Grid
    points::Vector{Point}

end

function PerturbedLatticeModel(h::AbstractHamiltonian, move::AbstractMoveModel, grid::Grid)
    pl = PerturbedLatticeModel(h, move, grid, deepcopy(grid.points))
    return pl
end

PerturbedLatticeModel(h::AbstractHamiltonian, move::AbstractMoveModel; N::Int=20, d::Int=2) = PerturbedLatticeModel(h, move, Grid(N, d))

Base.length(pl::PerturbedLatticeModel) = length(pl.grid)

function Base.show(io::IO, pl::PerturbedLatticeModel)
    print(io, "PerturbedLatticeModel  (d=$(pl.grid.d), N=$(pl.grid.N), $(length(pl)) points)")
end
