"""
Energy calculation functions for point processes.
"""


function compute_adjacency_matrix(points::Vector{Vector{Float64}}, R::Float64, d::Int)
    n_points = length(points)
    adjacency = zeros(Int, n_points, n_points)

    for i in 1:n_points
        for j in (i+1):n_points
            if d == 2
                dx = points[i][1] - points[j][1]
                dy = points[i][2] - points[j][2]
                dist_sq = dx*dx + dy*dy
            elseif d == 3
                dx = points[i][1] - points[j][1]
                dy = points[i][2] - points[j][2]
                dz = points[i][3] - points[j][3]
                dist_sq = dx*dx + dy*dy + dz*dz
            end

            if dist_sq <= R^2
                adjacency[i, j] = 1
                adjacency[j, i] = 1
            end
        end
    end

    return adjacency
end


function refresh_adjacency_matrix(adjacency::Matrix{Int}, points::Vector{Vector{Float64}}, R::Float64, new_point::Vector{Float64}, i::Int, d::Int)
    n_points = length(points)

    # Réinitialiser la ligne et la colonne i
    adjacency[i, :] .= 0
    adjacency[:, i] .= 0

    # Recalculer les adjacences pour le point i à sa nouvelle position
    for j in 1:n_points
        if j != i
            if d == 2
                dx = new_point[1] - points[j][1]
                dy = new_point[2] - points[j][2]
                dist_sq = dx*dx + dy*dy
            elseif d == 3
                dx = new_point[1] - points[j][1]
                dy = new_point[2] - points[j][2]
                dz = new_point[3] - points[j][3]
                dist_sq = dx*dx + dy*dy + dz*dz
            end
            if dist_sq <= R^2
                adjacency[i, j] = 1
                adjacency[j, i] = 1
            end
        end
    end

    return adjacency
end

function local_energy(adjacencyS::Matrix{Int}, adjacencyHC::Matrix{Int}, i::Int)
    if any(adjacencyHC[i, :] .== 1)
        return Inf
    end
    return Float64(sum(adjacencyS[i, :]))
end









"""
    local_energy_MC(pl::PerturbedLatticeModel)

Compute Monte Carlo estimates of local energy for all points.

This function fills `pl.loc_en_mc` with energy estimates for each point
across multiple perturbation samples.

"""
function local_energy_MC(pl::PerturbedLatticeModel)
    for i in 1:(2*pl.N+1)^pl.d
        for j in 1:pl.NMC
            new_config = deepcopy(pl.points)
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
                x = [rand(pl.rng, Uniform(pl.BorneUnif[k, 1], pl.BorneUnif[k, 2]))
                     for k in 1:pl.d]
                new_config[i] = pl.grid[i] .+ x
            end
            pl.loc_en_mc[i,j] = local_energy(pl, new_config, i)
        end
    end
    return nothing
end

"""
    estimate_int_DLR(pl::PerturbedLatticeModel, i::Int, beta::Float64)

Estimate the integral term in the DLR equation using Monte Carlo.
"""

function estimate_int_DLR(pl::PerturbedLatticeModel, i::Int, beta::Float64)
    num = 0.0
    denom = 0.0
    for j in 1:pl.NMC
        loc_en = pl.loc_en_mc[i,j]
        if loc_en == Inf
            num += 0
            denom += 0
        else
            num += loc_en*exp(-beta*loc_en)
            denom += exp(-beta*loc_en)
        end
    end
    if denom != 0
        return num / denom
    end
    return 0.0
end