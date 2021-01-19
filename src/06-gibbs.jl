"""
    Update the clusters using a Rao-Blackwellized Gibbs algorithm,
    see Sudderth (2006, p. 94, algorithm 2.2)
"""
function update_z!(c::ChainState, s::SuffStats{A, B, C}) where {A, B, C}
    @unpack refresh_rate, α, z, τ, O, U, rng = c
    @unpack n = s
    randperm!(rng, τ)

    counter = 0
    @inbounds @fastmath for i ∈ τ
        # Recompute from scratch each 30 iterations
        (counter += 1) % refresh_rate == 0 && suffstats0!(s, c)

        # Initialize the search
        z̄ = z[i]
        p̄ = -Inf
        k̄ = first(U)
        
        # Update zi
        for k ∈ O
            p = log_pl(s, c, i, k) + log(n[k] - (z̄ == k))
            p = p - log(-log(rand(rng)))
            p < p̄ && continue
            z[i] = k
            p̄ = p
        end
        p = log_pl(s, c, i, k̄) + log(α[1])
        p = p - log(-log(rand(rng)))
        p > p̄ && (z[i] = k̄)

        # Update the sufficient statistics
        z̄ != z[i] && suffstats!(s, c, z̄, z[i], i)
    end
end

function update_θ!(c::ChainState, s::SuffStats{A, B, C}) where {A, B, C}
    @unpack J, u, ν, r, S, Σ, μ = s
    @unpack O, rng = c
    for k in O, j ∈ 1:J
        Σ[j, k] = rand(rng, InverseWishart(ν[j, k], S[j, k]))
        Σ[j, k] = (Σ[j, k] + Σ[j, k]') / 2
        μ[j, k] = rand(rng, MvNormal(u[j, k], Σ[j, k] / r[j, k]))
    end
end

function update_γ!(
    c::ChainState, 
    s::SuffStats{A, B, C}, 
    pγ0::OffsetArray{Float64, 1, Vector{Float64}}
) where {A, B, C}
    @unpack n, J = s
    @unpack γ, O, rng = c

    # Resample γ[g], given the other γ's
    for g = 2:J
        # log-odds (numerator)
        γ[g] = 1
        suffstats0!(s, c)
        log_num = log(pγ0[sum(γ)])
        for k ∈ O, j ∈ (1, g)
            log_num += log_ml(s, j, k)
        end

        # log-odds (denominator)
        γ[g] = 0
        suffstats0!(s, c)
        log_den = log(pγ0[sum(γ)])
        for k ∈ O, j ∈ (1)
            log_den += log_ml(s, j, k)
        end

        # log-odds and new γ[g]
        log_odds = log_num - log_den
        γ[g] = rand(rng) <= 1 / (1 + exp(-log_odds))
    end
end

function update_α!(c::ChainState, s::SuffStats{A, B, C}) where {A, B, C}
    @unpack m, N = s
    @unpack a0, b0 = m
    @unpack α, K, rng = c

    ϕ = rand(rng, Beta(α[1] + 1.0, N))
    ψ = (a0 + K[1] - 1.0) / (N * (b0 - log(ϕ))); ψ = ψ / (1 + ψ)
    α[1] = rand(rng, Gamma(a0 + K[1] - (rand(rng) > ψ), 1.0 / (b0 - log(ϕ))))
end

function fit(
    m::Model{A, B},
    y::Vector{C},
    x::Vector{Int};
    K::Int = 5,
    iter::Int = 4000, 
    warmup::Int = iter ÷ 2,
    refresh_rate::Int = length(x),
    rng::MersenneTwister = MersenneTwister()
) where {A, B, C}
    # Initialization
    s = SuffStats(m = m, y = y, x = x)
    c = ChainState(N = s.N, J = s.J, rng = rng, K = [K])
    pγ1 = zeros(s.J - 1)
    pγ0 = ph0(s.J - 1, 1.0)

    # Updating
    @fastmath for t ∈ 1:iter
        suffstats0!(s, c)
        update_z!(c, s)
        update_α!(c, s)
        update_γ!(c, s, pγ0)
        t > warmup || continue
        pγ1[γcode(c.γ)] += 1
    end
    return pγ1 / (iter - warmup)
end

function fit(
    m::Model{A, B},
    y::Vector{C},
    x::Vector{Int},
    grid::Vector{Float64};
    K::Int = 5,
    iter::Int = 4000, 
    warmup::Int = iter ÷ 2,
    rng::MersenneTwister = MersenneTwister()
) where {A, B, C}
    # Initialization
    D = length(y[1])
    J = length(unique(x))
    s = SuffStats(m = m, y = y, x = x)
    c = ChainState(N = s.N, J = s.J, rng = rng, K = [K])
    pγ1 = zeros(s.J - 1)
    pγ0 = ph0(s.J - 1, 1.0)
    ygrid = Iterators.product(fill(grid, 2)...)
    ygrid = collect.(ygrid)[:] |> x -> hcat(x...)
    fgrid = [zeros(size(ygrid, 2)) for j = 1:J, z1 in 1:D, z2 in 1:D]
    M = size(ygrid, 2)

    # Updating
    @unpack γ, α, O = c
    @unpack n, N, μ, Σ = s
    @unpack r0, ν0 = m
    r1 = r0 + 1
    ν1 = ν0 + 1
    u0 = zeros(2)
    u1 = zeros(2)
    S0 = Cholesky(m.S0.factors[1:2, 1:2], :U, 0)
    S1 = deepcopy(S0)
    ggrid = zeros(M)
    yi = zeros(2)
    @fastmath for t ∈ 1:iter
        suffstats0!(s, c)
        update_z!(c, s)
        update_θ!(c, s)
        update_α!(c, s)
        update_γ!(c, s, pγ0)
        t > warmup || continue
        pγ1[γcode(c.γ)] += 1

        # Compute a base contribution to f(y)
        for i in 1:M
            S1.factors .= S0.factors
            yi .= ygrid[:, i]
            u1 .= (r0 * u0 + yi) / r1
            u1 .= (yi - u1) * √(r1 / r0)
            lowrankupdate!(S1, u1)
            ggrid[i] =
                0.5ν0 * logdet(S0) -
                0.5ν1 * logdet(S1) +
                logmvgamma(2, ν1 / 2) -
                logmvgamma(2, ν0 / 2) +
                0.5 * 2 * log(r0 / r1) -
                0.5 * 2 * log(π) * (r1 - r0)
            ggrid[i] = exp(ggrid[i])
        end

        # Accumulate f(y)
        μsub = [0.0; 0.0]
        Σsub = [1.0 0.0; 0.0 1.0]
        for j = 1:J, z1 = 1:D, z2 = z1+1:D
            for k in O
                μsub .= μ[j^γ[j], k][[z1; z2]]
                Σsub .= Σ[j^γ[j], k][[z1; z2], [z1; z2]]
                fgrid[j, z1, z2] += pdf(MvNormal(μsub, Σsub), ygrid) * n[k] / (N + α[1])
            end
            fgrid[j, z1, z2] += ggrid * α[1] / (N + α[1])
        end
    end

    # Convert fgrid into a dataframe
    dfs = [
        DataFrame(
            j    = j, 
            var1 = z1, 
            var2 = z2, 
            y1   = ygrid[1, :],
            y2   = ygrid[2, :],
            f    = fgrid[j, z1, z2] / (iter - warmup)
        ) for j = 1:J, z1 = 1:D, z2 = 1:D
    ]
    df = reduce(vcat, dfs)
    filter!([:var1, :var2] => (x, y) -> x < y, df)
    return pγ1 / (iter - warmup), df
end

function train(
    y::Matrix{Float64}, 
    x::Vector{Int}; 
    r0::Int = 1, 
    v0::Int = size(y, 2) + 2, 
    u0::Vector{Float64} = zeros(size(y, 2)), 
    S0::Matrix{Float64} = Matrix{Float64}(I(size(y, 2))), 
    a0::Float64 = 1.0, 
    b0::Float64 = 1.0, 
    z0::Float64 = 1.0,
    iter::Int   = 4000,
    warmup::Int = iter ÷ 2,
    rng::MersenneTwister = MersenneTwister()
)
    N, D = size(y)
    J = length(unique(x))
    m = MANOVABNPTest.Model(
        D  = D,
        r0 = r0,
        ν0 = v0,
        u0 = u0,
        S0 = cholesky(S0),
        a0 = a0,
        b0 = b0,
        ζ0 = z0
    )
    y = standardize(ZScoreTransform, y, dims = 1)
    y = [SVector{D}(y[i, :]) for i ∈ 1:N]
    ps = MANOVABNPTest.fit(m, y, x; iter = iter, warmup = warmup, rng = rng)
    γs = [γvector(J, u)[2:end] for u in 1:length(ps)]
    DataFrame(hypothesis = γs, prob = ps)
end