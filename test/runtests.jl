using MANOVABNPTest
using StatsBase, LinearAlgebra, StaticArrays

function test_sample(h::Int, l::Int, H0::Int)
    N = 4 * [50, 150, 300]
    K = [
        0.625 1.498 1.162;
        0.819 1.806 1.359;
        1.051 1.994 1.652
    ]
    γ = γvector(4, H0)
    Σ = cholesky([1.0 0.3; 0.3 1.0])
    x = 1 .+ (0:N[h]-1) .% 4
    y = rand(N[h], 2) * Σ.U
    for i in 1:N[h]
        x[i] == 2 && (y[i, :] .+= K[l, 1] * γ[2])
        x[i] == 3 && (y[i, :] .*= K[l, 2] ^ γ[3])
        x[i] == 4 && (y[i, :] .+= K[l, 3] * γ[4] * rand([-1, 1]))
    end
    ȳ, S = mean_and_cov(y)
    y = (y .- ȳ) / cholesky(S).U
    return y, x
end
y, x = test_sample(1, 1, 8);
println(train(y, x))
1+1