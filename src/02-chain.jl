"""
    Chain state in our Gibbs sampler
"""
@with_kw struct ChainState
    N::Int
    J::Int
    K̄::Int = 64
    K::Vector{Int}       = [5]
    rng::MersenneTwister = Random.MersenneTwister()
    τ::Vector{Int}       = randperm(rng, N)
    γ::Vector{Int}       = [0; ones(Int, J - 1)]
    z::Vector{Int}       = [1:K[1]; rand(rng, 1:K[1], N - K[1])]
    U::Set{Int}          = Set(K[1]+1:K̄-1)
    O::Set{Int}          = Set(1:K[1])
    α::Vector{Float64}   = [1.0]
end