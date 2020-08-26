module GGA_MANOVA_DDP

# Imports
using Distributions: Beta, Gamma
using LinearAlgebra: cholesky, I, lowrankupdate!, lowrankdowndate!, logdet
using OffsetArrays: OffsetArray
using Parameters: @with_kw, @unpack
using StaticArrays: MMatrix, MVector, @MVector, @MMatrix
using StatsFuns: logmvgamma
using Random: shuffle!

# Write your package code here.
include("utils.jl")
include("chain.jl")
include("model.jl")
include("suffstats.jl")
include("logl.jl")
include("gibbs.jl")

# Exports
export Model, SuffStats, suffstats!, ChainState, log_pl, log_ml, ph0, fit
export update_z!, update_γ!, update_α!
export γcode, γvector

end
