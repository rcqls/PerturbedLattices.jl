struct Grid <: AbstractLattice
    N::Int
    d::Int
    points::Vector{Point}

    
    function Grid(N::Int, d::Int) 
        gr = new(N, d, Vector{Point}(undef, (2*N+1)^d))
        init!(gr)
        return gr
    end
end

function Base.length(g::Grid)
    return (2*g.N+1)^g.d
end

function Base.size(g::Grid)
    return tuple(repeat([2*g.N+1], g.d))
end

function init!(g::Grid)
    l = length(g)

    if g.d == 2
        for i in 0:(2*g.N)
            for j in 0:(2*g.N)
                g.points[(i + 1) + j * (2 * g.N + 1)] = [Float64(i - g.N), Float64(j - g.N)]
            end
        end
    elseif g.d == 3
        for i in 0:(2*g.N)
            for j in 0:(2*g.N)
                for k in 0:(2*g.N)
                    g.points[(i + 1) + j * (2 * g.N + 1) + k * (2 * g.N + 1)^2] = [Float64(i - g.N), Float64(j - g.N), Float64(k - g.N)]
                end
            end
        end
    else
        throw(ArgumentError("Dimension d must be 2 or 3"))
    end

end


#SubGrid is just a grid with SubGrid.N < Grid.N

