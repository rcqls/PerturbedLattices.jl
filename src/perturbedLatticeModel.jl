"""
    PerturbedLatticeModel

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

mutable struct PerturbedLattice
    N::Int
    d::Int
    NMC::Int

    H::Union{String, Nothing}
    RS::Union{Float64, Nothing}
    RHC::Union{Float64, Nothing}
    beta::Union{Float64, Nothing}

    Distrib::Union{String, Nothing}
    BorneUnif::Union{Matrix{Float64}, Nothing}
    Cov::Union{Matrix{Float64}, Nothing}

    grid::Vector{Vector{Float64}}
    points::Vector{Vector{Float64}}
    adjacencyS::Matrix{Int}
    adjacencyHC::Matrix{Int}

    seed::Union{Int, Nothing}
    rng::AbstractRNG

    function PerturbedLattice(d::Int, N::Int;
                              H::Union{String, Nothing}=nothing,
                              RS::Union{Real, Nothing}=nothing,
                              RHC::Union{Real, Nothing}=nothing,
                              beta::Union{Real, Nothing}=nothing,
                              sigma::Union{Real, Nothing}=nothing,
                              Distrib::Union{String, Nothing}=nothing,
                              BorneUnif::Union{Matrix, Nothing}=nothing,
                              seed::Union{Int, Nothing}=nothing,
                              NMC::Int64=10000)

        # Type conversions
        RS_float = isnothing(RS) ? nothing : Float64(RS)
        RHC_float = isnothing(RHC) ? nothing : Float64(RHC)
        beta_float = isnothing(beta) ? nothing : Float64(beta)
        sigma_float = isnothing(sigma) ? nothing : Float64(sigma)
        BorneUnif_mat = isnothing(BorneUnif) ? nothing : Float64.(BorneUnif)
        Cov_mat = isnothing(Cov) ? nothing : Float64.(Cov)

        # Initialize random number generator
        rng = isnothing(seed) ? Random.default_rng() : MersenneTwister(seed)

        # Create instance without grid, points, and loc_en
        obj = new(N, d, NMC, H, RS_float, RHC_float, beta_float,sigma_float,
                  Distrib, BorneUnif_mat, Vector{Vector{Float64}}(), Vector{Vector{Float64}}(),
                  Float64[], zeros(0,0), seed, rng)

        # Initialize grid and points
        obj.grid = create_grid(obj)
        obj.points = deepcopy(obj.grid)
        obj.loc_en = local_energy_vec(obj)
        obj.adjacencyS = compute_adjacency_matrix(obj.points, obj.RS)
        obj.adjacencyHC = compute_adjacency_matrix(obj.points, obj.RHC)

        return obj
    end
end

# Extend Base functions for convenience
Base.length(pl::PerturbedLattice) = length(pl.points)

function Base.show(io::IO, pl::PerturbedLattice)
    print(io, "PerturbedLattice(d=$(pl.d), N=$(pl.N), $(length(pl)) points)")
end