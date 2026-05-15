"""
Monte Carlo simulation functions using Metropolis-Hasting algorithm. 
"""
function iterate!(pl::PerturbedLatticeV1)
    new_config = pl.points
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
    
    old_loc_en = local_energy(pl.adjacencyS, pl.adjacencyHC, i)
    new_adjacencyS = refresh_adjacency_matrix(copy(pl.adjacencyS), pl.points, pl.RS, new_config[i], i, pl.d)
    new_adjacencyHC = isnothing(pl.RHC) ? copy(pl.adjacencyHC) : refresh_adjacency_matrix(copy(pl.adjacencyHC), pl.points, pl.RHC, new_config[i], i, pl.d)
    new_loc_en = local_energy(new_adjacencyS, new_adjacencyHC, i)

    # Metropolis-Hastings acceptance ratio
    if new_loc_en == Inf && old_loc_en == Inf
        r = 1.0
    elseif new_loc_en == Inf
        r = 0.0
    elseif old_loc_en == Inf
        r = 1.0
    else
        r = exp(-pl.beta * (new_loc_en - old_loc_en))
    end

    if rand(pl.rng) <= r
        pl.points = new_config
        pl.adjacencyS = new_adjacencyS
        pl.adjacencyHC = new_adjacencyHC
    end
    return nothing
end
