using DataFrames
using GGA_MANOVA_DDP
using LinearAlgebra
using StaticArrays
using StatsBase
using RCall
using JLD2

function γcode(γ)
    code = 1
    @inbounds for j = 2:length(γ)
        code += γ[j] * 2^(j - 2)
    end
    code
end

function TestParameters(h::Int, l::Int, c::Int)
    # Parameters for the reference group
    μ0 = zeros(2)
    Σ0 = [1.0 0.3; 0.3 1.0]

    # Parameters for the other groups
    c1 = (0.625, 0.819, 1.051)
    c2 = (1.498, 1.806, 1.994)
    c3 = (1.162, 1.359, 1.652)
    μ = [[copy(μ0) for j ∈ 1:4, k ∈ 1:2] for l ∈ 1:3, c ∈ 1:8]
    Σ = [[copy(Σ0) for j ∈ 1:4, k ∈ 1:2] for l ∈ 1:3, c ∈ 1:8]
    for l ∈ 1:3, c ∈ 1:8, k ∈ 1:2
        γ1, γ2, γ3, γ4 = γvector(4, c)
        μ[l, c][2, k] .+= c1[l] * γ2
        Σ[l, c][3, k] .*= c2[l] ^ γ3
        μ[l, c][4, k] .+= c3[l] * γ4 * 2 * (k - 1.5)
    end

    # Sample sizes
    N = [200, 600, 900]
    return μ[l, c], cholesky.(Σ[l, c]), N[h]
end

R"""
library(car)
pillai_pval <- function(x) {
	# the following code is adapted from summary.manova
	eigs <- Re(eigen(qr.coef(qr(x$SSPE), x$SSPH), symmetric = FALSE)$values)
	stat <- stats:::Pillai(eigs, x$df, x$df.residual)
	pval <- pf(stat[2], stat[3], stat[4], lower.tail = FALSE)
	return(pval)
}
"""

R"""
get_hypothesis <- function(df) {
    fit <- lm(cbind(y1, y2) ~ x + 0, df)
    pvals <- rep(0, 3)
    for (j in 1:3) {
        x <- linearHypothesis(fit, h = matrix(c(-1, rep(1:3) == j), 1, 4))
        pvals[j] <- pillai_pval(x)
    }
    gvector <- c(0, pvals <= 0.05 / 3)
    return(gvector)
}
"""

function test_sample(h::Int, l::Int, c::Int)
    μ, Σ, N = TestParameters(h, l, c)
    x = [1 + (i - 1) % 4 for i ∈ 1:N]
    z = [1 + (rand() <= 0.5) for i ∈ 1:N]
    y = [μ[x[i], z[i]] + Σ[x[i], z[i]].L * randn(2) for i ∈ 1:N]
    y = Vector{Float64}.(y)
    y = standardize(ZScoreTransform, hcat(y...)', dims = 1)
    df = DataFrame(y1 = y[:, 1], y2 = y[:, 2], x = CategoricalArray(x))
    return df
end

function mc_experiment_02()
    nsuccess = zeros(3, 3, 8, 8)
    for h ∈ 1:3, l ∈ 1:3, c ∈ 1:8
        for sim ∈ 1:100
            df = test_sample(h, l, c);        @rput df
            R"gvector <- get_hypothesis(df)"; @rget gvector
            nsuccess[h, l, c, Int(γcode(gvector))] += 1
        end
    end
    return nsuccess
end

tbl = mc_experiment_02()
@save "data/mc_experiment_02.jld2" tbl