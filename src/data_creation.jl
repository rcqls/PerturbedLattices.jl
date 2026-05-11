""" 
    Extraction of points within a window and separation between interior and boundary points. 
"""

function points_in_window(pl::PerturbedLattice, W::Matrix{Float64}, m_n::Float64)
    points_D_n = Vector{Vector{Float64}}()
    grid_D_n = Vector{Vector{Float64}}()
    points_W_n_boundary = Vector{Vector{Float64}}()

    
    if pl.d ==2
        for i in 1:(2*pl.N+1)^2
            x, y = pl.points[i][1], pl.points[i][2]
            if W[1,1] + pl.RS <= x <= W[1,2] - pl.RS && W[2,1] + pl.RS <= y <= W[2,2] - pl.RS
                if W[1,1] + pl.RS + m_n  <= pl.grid[i][1] <= W[1,2] - pl.RS - m_n  && W[2,1] + pl.RS + m_n <= pl.grid[i][2] <= W[2,2] - pl.RS - m_n
                    push!(points_D_n, [x, y])
                    push!(grid_D_n, [pl.grid[i][1], pl.grid[i][2]])
                else 
                    push!(points_W_n_boundary, [x, y])
                end
            elseif W[1,1] <= x <= W[1,2] && W[2,1] <= y <= W[2,2]
                push!(points_W_n_boundary, [x, y])
            end
        end
    end

    return points_D_n, grid_D_n, points_W_n_boundary
end