mutable struct Points 
    grid::Grid
    points::Vector{Point}
end

function Points(N::Int, d::Int) #N et d définissent la fenêtre d'observation
    g = Grid(N, d)
    points = Vector{Point}(undef, (2*N+1)^d)
    return Points(g, points)
end