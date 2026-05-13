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

    function PerturbedLattice(N::Int;
                              d::Int=2,
                              H::Union{String, Nothing}=nothing,
                              RS::Union{Real, Nothing}=nothing,
                              RHC::Union{Real, Nothing}=nothing,
                              beta::Union{Real, Nothing}=nothing,
                              sigma::Union{Real, Nothing}=nothing,
                              Distrib::Union{String, Nothing}=nothing,
                              BorneUnif::Union{Matrix, Nothing}=nothing,
                              seed::Union{Int, Nothing}=nothing,
                              NMC::Int64=10000)

        RS_float = isnothing(RS) ? nothing : Float64(RS)
        RHC_float = isnothing(RHC) ? nothing : Float64(RHC)
        beta_float = isnothing(beta) ? nothing : Float64(beta)
        BorneUnif_mat = isnothing(BorneUnif) ? nothing : Float64.(BorneUnif)

        if !isnothing(sigma)
            sigma_f = Float64(sigma)
            Cov_mat = sigma_f^2 * Matrix{Float64}(I, d, d)
            Distrib = isnothing(Distrib) ? "Gauss" : Distrib
        else
            Cov_mat = nothing
        end

        rng = isnothing(seed) ? Random.default_rng() : MersenneTwister(seed)

        obj = new(N, d, NMC, H, RS_float, RHC_float, beta_float,
                  Distrib, BorneUnif_mat, Cov_mat,
                  Vector{Vector{Float64}}(), Vector{Vector{Float64}}(),
                  zeros(Int, 0, 0), zeros(Int, 0, 0),
                  seed, rng)

        obj.grid = create_grid(obj)
        obj.points = deepcopy(obj.grid)
        n = length(obj.points)
        obj.adjacencyS = isnothing(RS_float) ? zeros(Int, n, n) : compute_adjacency_matrix(obj.points, RS_float, d)
        obj.adjacencyHC = isnothing(RHC_float) ? zeros(Int, n, n) : compute_adjacency_matrix(obj.points, RHC_float, d)

        return obj
    end
end

Base.length(pl::PerturbedLattice) = length(pl.points)

function Base.show(io::IO, pl::PerturbedLattice)
    print(io, "PerturbedLattice(d=$(pl.d), N=$(pl.N), $(length(pl)) points)")
end
