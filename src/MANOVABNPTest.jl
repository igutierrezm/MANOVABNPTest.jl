module MANOVABNPTest

# Imports
using DataFrames: DataFrame, filter!
using Distributions: Beta, Gamma, InverseWishart, MvNormal, pdf
using LinearAlgebra: cholesky, Cholesky, I, lowrankupdate!, lowrankdowndate!, logdet, Symmetric
using OffsetArrays: OffsetArray
using Parameters: @with_kw, @unpack
using StaticArrays: MMatrix, MVector, SVector, @MVector, @MMatrix
using StatsFuns: logmvgamma
using StatsBase: mean, cov, standardize, ZScoreTransform
using Random: randperm, randperm!, MersenneTwister

# Write your package code here.
include("01-utils.jl")
include("02-chain.jl")
include("03-model.jl")
include("04-suffstats.jl")
include("05-logl.jl")
include("06-gibbs.jl")

# Exports
export Model, SuffStats, suffstats0!, suffstats!, ChainState, log_pl, log_ml, ph0, fit, train, train_γ
export update_z!, update_γ!, update_α!
export γcode, γvector, γstr, γvec
export test_sample

end
