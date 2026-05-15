abstract type AbstractHamiltonian end

mutable struct StraussHamiltonian <: AbstractHamiltonian
    beta::Float64
    radius::Float64 
    # fields required for AbstractHamiltonian
    lattice::AbstractLattice
    adjacency::Matrix{Int}

    function StraussHamiltonian(beta::Float64, radius::Float64, lattice::AbstractLattice)
        sh = new(beta, radius, lattice)
        adjacency_matrix!(sh)
        return sh
    end
end


mutable struct HardCoreHamiltonian <: AbstractHamiltonian
    radius::Float64
    # fields required for AbstractHamiltonian
    lattice::AbstractLattice
    adjacency::Matrix{Int}

    function HardCoreHamiltonian(radius::Float64, lattice::AbstractLattice)
        hch = new(radius, lattice)
        adjacency_matrix!(hch)
        return hch
    end

end


function local_energy(h::StraussHamiltonian, i::Int)
    return Float64(h.beta*sum(h.adjacency[i, :]))
end

function local_energy(h::HardCoreHamiltonian, i::Int)
    if any(h.adjacency[i, :] .== 1)
        return Inf
    end
    return 0.0
end

function adjacency_matrix!(h::AbstractHamiltonian)
    points = h.lattice.points
    d = h.lattice.d
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

            if dist_sq <= h.radius^2
                adjacency[i, j] = 1
                adjacency[j, i] = 1
            end
        end
    end

    h.adjacency = adjacency
end


function adjacency_matrix!(h::AbstractHamiltonian, i::Int, new_point::Point)
    n_points = length(h.lattice.points)
    d = h.lattice.d
    points = h.lattice.points
    adjacency = h.adjacency

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
            if dist_sq <= h.radius^2
                adjacency[i, j] = 1
                adjacency[j, i] = 1
            end
        end
    end

end