using DataFrames, MANOVABNPTest, StatsBase, LinearAlgebra, StaticArrays, Random
rng = MersenneTwister(1)

function test_sample_01(D::Int, h::Int, l::Int, H0::Int)
    N = 4 * [50, 150, 300]
    K = [
        0.625 1.498 1.162;
        0.819 1.806 1.359;
        1.051 1.994 1.652
    ]
    γ = γvector(4, H0)
    Σ = cholesky(collect(I(D)))
    x = 1 .+ (0:N[h]-1) .% 4
    y = rand(N[h], D) * Σ.U
    for i in 1:N[h]
        x[i] == 2 && (y[i, :] .+= K[l, 1] * γ[2])
        x[i] == 3 && (y[i, :] .*= K[l, 2] ^ γ[3])
        x[i] == 4 && (y[i, :] .+= K[l, 3] * γ[4] * rand([-1, 1]))
    end
    ȳ, S = mean_and_cov(y)
    y = (y .- ȳ) / cholesky(S).U
    y = [SVector{D}(y[i, :]) for i ∈ 1:N[h]]
    return y, x
end

function test_sample_02(D::Int, h::Int, l::Int, H0::Int)
    N = 4 * [50, 150, 300]
    K = [
        0.625 1.498 1.162;
        0.819 1.806 1.359;
        1.051 1.994 1.652
    ]
    γ = γvector(4, H0)
    Σ = cholesky(collect(I(D)))
    x = 1 .+ (0:N[h]-1) .% 4
    y = rand(N[h], D) * Σ.U
    for i in 1:N[h]
        x[i] == 2 && (y[i, :] .+= K[l, 1] * γ[2])
        x[i] == 3 && (y[i, :] .*= K[l, 2] ^ γ[3])
        x[i] == 4 && (y[i, :] .+= K[l, 3] * γ[4] * rand([-1, 1]))
    end
    return y, x
end

# Exp 1
D = 2
m = MANOVABNPTest.Model(D = D)
y, x = test_sample_01(D, 1, 1, 2)
pγ = MANOVABNPTest.fit(m, y, x; iter = 4000, refresh_rate = 50, rng = rng)
println(pγ)

# Exp 2
D = 2
m = MANOVABNPTest.Model(D = D)
y, x = test_sample_01(D, 1, 1, 3)
pγ = MANOVABNPTest.fit(m, y, x; iter = 4000, refresh_rate = 50, rng = rng)
println(pγ)

# Exp 3
D = 10
m = MANOVABNPTest.Model(D = D)
y, x = test_sample_01(D, 1, 1, 2)
pγ = MANOVABNPTest.fit(m, y, x; iter = 4000, refresh_rate = 50, rng = rng)
println(pγ)

# Exp 4
D = 10
m = MANOVABNPTest.Model(D = D)
y, x = test_sample_01(D, 1, 1, 3)
pγ = MANOVABNPTest.fit(m, y, x; iter = 4000, refresh_rate = 50, rng = rng)
println(pγ)

# Exp 5
D = 10
m = MANOVABNPTest.Model(D = D)
y, x = test_sample_01(D, 1, 2, 2)
pγ = MANOVABNPTest.fit(m, y, x; iter = 4000, refresh_rate = 50, rng = rng)
println(pγ)

# Exp 6
D = 10
m = MANOVABNPTest.Model(D = D)
y, x = test_sample_01(D, 1, 2, 3)
pγ = MANOVABNPTest.fit(m, y, x; iter = 4000, refresh_rate = 50, rng = rng)
println(pγ)

# Exp 7
D = 10
m = MANOVABNPTest.Model(D = D)
y, x = test_sample_01(D, 1, 3, 2)
pγ = MANOVABNPTest.fit(m, y, x; iter = 4000, refresh_rate = 50, rng = rng)
println(pγ)

# Exp 8
D = 10
m = MANOVABNPTest.Model(D = D)
y, x = test_sample_01(D, 1, 3, 3)
pγ = MANOVABNPTest.fit(m, y, x; iter = 4000, refresh_rate = 50, rng = rng)
println(pγ)

# grid = LinRange(-3, 3, 10) |> collect
# m = MANOVABNPTest.Model(D = 2)
# MANOVABNPTest.fit(m, y, x, grid; iter = 200, rng = rng);

# Exp 9
y, x = test_sample_01(2, 1, 1, 2)
Y = Matrix{Float64}(vcat(y'...))
println(train(Y, x))
