"""
Grid creation and manipulation functions.
"""


function create_grid(pl::PerturbedLatticeModel)
    grid = Vector{Vector{Float64}}()

    if pl.d == 2
        for i in 0:(2*pl.N)
            for j in 0:(2*pl.N)
                push!(grid, [Float64(i - pl.N), Float64(j - pl.N)])
            end
        end
    elseif pl.d == 3
        for i in 0:(2*pl.N)
            for j in 0:(2*pl.N)
                for k in 0:(2*pl.N)
                    push!(grid, [Float64(i - pl.N), Float64(j - pl.N), Float64(k - pl.N)])
                end
            end
        end
    else
        throw(ArgumentError("Dimension d must be 2 or 3"))
    end

    return grid
end

"""
    shift!(pl::PerturbedLatticeModel, u::Vector)

Translate all points by vector u.

"""
function shift!(pl::PerturbedLatticeModel, u::Vector)
    for i in 1:length(pl.points)
        pl.points[i] .+= u
    end
    return nothing
end
