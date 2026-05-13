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

function local_energy(adjacency::Matrix{Int}, i::Int)
    return Float64(sum(adjacency[i, :]))
end









