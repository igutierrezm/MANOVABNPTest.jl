using Distributed, Test, JLD2
nworkers() == 8 || addprocs(8)
nworkers()

@everywhere begin
    using GGA_MANOVA_DDP
    using LinearAlgebra
    using StaticArrays
    using StatsBase
end

@everywhere function TestParameters(h::Int, l::Int, c::Int)
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

@everywhere function test_sample(h::Int, l::Int, c::Int)
    μ, Σ, N = TestParameters(h, l, c)
    x = [1 + (i - 1) % 4 for i ∈ 1:N]
    z = [1 + (rand() <= 0.5) for i ∈ 1:N]
    y = [μ[x[i], z[i]] + Σ[x[i], z[i]].L * randn(2) for i ∈ 1:N]
    y = Vector{Float64}.(y)
    y = standardize(ZScoreTransform, hcat(y...)', dims = 1)
    y = SVector{2, Float64}.([y[i, :] for i ∈ 1:N])
    return y, x
end

@everywhere function experiment(h::Int, l::Int, c::Int)
    y, x = test_sample(h, l, c)
    return GGA_MANOVA_DDP.fit(Model(D = 2), y, x; iter = 4000)
end

for h ∈ 2:3, l ∈ 1:3, c ∈ 1:8
    tbl = pmap(i -> experiment(h, l, c), 1:100);
    @save "../data/tbl_$(h)_$(l)_$(c).jld2" tbl
end

@testset "GGA_MANOVA_DDP.jl" begin
    # # Write your tests here.

    # N = 100
    # J = 4
    # D = 2

    # y = [randn(D)  for i = 1:N]
    # x = [1 + i % J for i = 1:N]
    # z = rand(1:2, N)
    # γ = ones(Int, J)

    # m = Model(D = 2)
    # c = ChainState(N = N, J = J)
    # s = SuffStats(m = m, y = y, x = x)
    # s = suffstats!(s, c)
    # # s = suffstats!(s, c, c.z[1], 6 - c.z[1], 1)
    # # log_pl(s, c, 1, 1)
    # # log_ml(s, 1, 1)
    # update_z!(c, s)

    # pγ0 = ph0(J - 1, 1.0)

    # N = 200
    # p = 0.5
    # μ = [@MVector [0.0, 0.0]         for j ∈ 1:4, k ∈ 1:2];
    # Σ = [@MMatrix [1.0 0.3; 0.3 1.0] for j ∈ 1:4, k ∈ 1:2] .|> cholesky;
    # for j ∈ 2:3, k ∈ 1:2
    #     @. μ[j, k] .+= 0.7
    # end
    # y, x, z = test_sample(μ, Σ, N, p);
    # m = Model(D = 2)
    # fit(m, y, x; iter = 4000)
end
