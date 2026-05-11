"""
Monte Carlo simulation functions.
"""

"""
    iterate!(pl::PerturbedLatticeModel)

Perform one iteration of the Metropolis-Hastings algorithm.

Randomly selects a point, proposes a new position based on the perturbation
distribution, and accepts or rejects based on the energy change.

"""
function iterate!(pl::PerturbedLatticeModel)
    new_config = deepcopy(pl.points)
    i = rand(pl.rng, 1:(2*pl.N+1)^pl.d)

    # Propose new position based on distribution
    if pl.Distrib == "Gauss"
        if pl.d == 2
            dist = MvNormal([0.0, 0.0], pl.Cov)
            x = rand(pl.rng, dist)
            new_config[i] = pl.grid[i] .+ x
        elseif pl.d == 3
            dist = MvNormal([0.0, 0.0, 0.0], pl.Cov)
            x = rand(pl.rng, dist)
            new_config[i] = pl.grid[i] .+ x
        end
    elseif pl.Distrib == "Unif"
        x = [rand(pl.rng, Uniform(pl.BorneUnif[j, 1], pl.BorneUnif[j, 2]))
             for j in 1:pl.d]
        new_config[i] = pl.grid[i] .+ x
    end
    
    #compute old and new local energy
    old_loc_en = local_energy(pl.adjacencyS, pl.adjacencyHC, i)
    new_adjacencyS = refresh_adjacency_matrix(pl.adjacencyS, pl.points, new_config[i], pl.RS, i, pl.d)
    new_adjacencyHC = refresh_adjacency_matrix(pl.adjacencyHC, pl.points, new_config[i], pl.RHC, i, pl.d)
    new_loc_en = local_energy(new_adjacencyS, new_adjacencyHC, i)

    # Metropolis-Hastings acceptance ratio
    if new_loc_en == Inf && old_loc_en == Inf
        r = 1.0
    elseif new_loc_en == Inf
        r = 0.0
    elseif old_loc_en == Inf
        r = 1.0 
    end
    r = exp(-pl.beta * (new_loc_en - old_loc_en))

    if rand(pl.rng) <= r
        pl.points = new_config
        pl.adjacencyS = new_adjacencyS
        pl.adjacencyHC = new_adjacencyHC
    end
    return nothing
end
