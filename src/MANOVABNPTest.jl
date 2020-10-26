module MANOVABNPTest

# Imports
using DataFrames: DataFrame, filter!
using Distributions: Beta, Gamma, InverseWishart, MvNormal, pdf
using LinearAlgebra: cholesky, I, lowrankupdate!, lowrankdowndate!, logdet
using OffsetArrays: OffsetArray
using Parameters: @with_kw, @unpack
using StaticArrays: MMatrix, MVector, SVector, @MVector, @MMatrix
using StatsFuns: logmvgamma
using Random: shuffle!

# Write your package code here.
include("01-utils.jl")
include("02-chain.jl")
include("03-model.jl")
include("04-suffstats.jl")
include("05-logl.jl")
include("06-gibbs.jl")

# Exports
export Model, SuffStats, suffstats!, ChainState, log_pl, log_ml, ph0, fit, train
export update_z!, update_γ!, update_α!
export γcode, γvector
export test_sample

end
