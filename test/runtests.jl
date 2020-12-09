using MANOVABNPTest, StatsBase, LinearAlgebra, StaticArrays, Random

function γvector(J, code)
    γ = zeros(Int, J)
    code -= 1
    for j = J:-1:2
        if (code ÷ 2^(j-2)) > 0
            code -= 2^(j-2)
            γ[j] = 1
        end
    end
    return γ
end

function test_sample_01(h::Int, l::Int, H0::Int)
    N = 4 * [50, 150, 300]
    K = [
        0.625 1.498 1.162;
        0.819 1.806 1.359;
        1.051 1.994 1.652
    ]
    γ = γvector(4, H0)
    Σ = cholesky(collect(I(2)))
    x = 1 .+ (0:N[h]-1) .% 4
    y = rand(N[h], 2) * Σ.U
    for i in 1:N[h]
        x[i] == 2 && (y[i, :] .+= K[l, 1] * γ[2])
        x[i] == 3 && (y[i, :] .*= K[l, 2] ^ γ[3])
        x[i] == 4 && (y[i, :] .+= K[l, 3] * γ[4] * rand([-1, 1]))
    end
    ȳ, S = mean_and_cov(y)
    y = (y .- ȳ) / cholesky(S).U
    y = [SVector{2}(y[i, :]) for i ∈ 1:N[h]]
    return y, x
end

function test_sample_02(h::Int, l::Int, H0::Int)
    N = 4 * [50, 150, 300]
    K = [
        0.625 1.498 1.162;
        0.819 1.806 1.359;
        1.051 1.994 1.652
    ]
    γ = γvector(4, H0)
    Σ = cholesky(collect(I(2)))
    x = 1 .+ (0:N[h]-1) .% 4
    y = rand(N[h], 2) * Σ.U
    for i in 1:N[h]
        x[i] == 2 && (y[i, :] .+= K[l, 1] * γ[2])
        x[i] == 3 && (y[i, :] .*= K[l, 2] ^ γ[3])
        x[i] == 4 && (y[i, :] .+= K[l, 3] * γ[4] * rand([-1, 1]))
    end
    return y, x
end

y, x = test_sample_01(1, 1, 1)
m = MANOVABNPTest.Model(D = 2)
rng = MersenneTwister(1)
ChainState(N = 10, J = 4, rng = rng)

pγ1 = MANOVABNPTest.fit(m, y, x; iter = 4000, rng = rng)
grid = LinRange(-3, 3, 10) |> collect
m = MANOVABNPTest.Model(D = 2)
MANOVABNPTest.fit(m, y, x, grid; iter = 200, rng = rng);

y, x = test_sample_02(1, 1, 6)
println(MANOVABNPTest.train(y, x, rng = rng))
