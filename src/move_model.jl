abstract type AbstractMoveModel end

struct GaussianMoveModel <: AbstractMoveModel
    cov::Matrix{Float64}
    d::Int64
end

struct UniformMoveModel <: AbstractMoveModel
    bounds::Vector{Vector{Float64}}
    d::Int64
end


function Base.rand(rng::AbstractRNG, model::GaussianMoveModel)
    return rand(rng, MvNormal(zeros(model.d), model.cov))
end

function Base.rand(rng::AbstractRNG, model::UniformMoveModel)
   return [rand(rng, Uniform(model.bounds[i][1], model.bounds[i][2])) for i in 1:model.d]
end