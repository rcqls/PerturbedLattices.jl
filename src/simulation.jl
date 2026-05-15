"""
Monte Carlo simulation functions using Metropolis-Hasting algorithm. 
"""
function iterate!(rng::AbstractRNG, pl::PerturbedLatticeModel)
    # Proposal of new point 
    i = rand(rng, 1:length(pl))
    new_point = pl.grid.points[i] .+ rand(rng, pl.move)
    
    old_adjacency = copy(pl.h.adjacency)
    
    old_loc_en = local_energy(pl.h, i)
    adjacency_matrix!(pl.h, i, new_point)
    new_loc_en = local_energy(pl.h, i)

    # Metropolis-Hastings acceptance ratio
    if new_loc_en == Inf && old_loc_en == Inf
        r = 1.0
    elseif new_loc_en == Inf
        r = 0.0
    elseif old_loc_en == Inf
        r = 1.0
    else
        r = exp(-(new_loc_en - old_loc_en))
    end

    if rand(rng) <= r
        pl.points[i] = new_point
    else 
        pl.h.adjacency = old_adjacency
    end
end

function Random.rand!(rng::AbstractRNG, pl::PerturbedLatticeModel; NMC::Int=1000)
    for i in 1:NMC
        iterate!(rng, pl)
    end
end

Random.rand!(pl::PerturbedLatticeModel; NMC::Int=1000) = rand!(Random.default_rng(), pl; NMC=NMC)