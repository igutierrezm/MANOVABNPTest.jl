"""
    Update the clusters using a Rao-Blackwellized Gibbs algorithm,
    see Sudderth (2006, p. 94, algorithm 2.2)
"""
function update_z!(c::ChainState, s::SuffStats{A, B, C}) where {A, B, C}
    @unpack α, z, τ, O, U = c
    @unpack n = s
    shuffle!(τ)

    @inbounds @fastmath for i ∈ τ
        # Initialize the search
        z̄ = z[i]
        p̄ = -Inf
        k̄ = first(U)
        
        # Update zi
        for k ∈ O
            p = log_pl(s, c, i, k) + log(n[k] - (z̄ == k))
            p = p - log(-log(rand()))
            p < p̄ && continue
            z[i] = k
            p̄ = p
        end
        p = log_pl(s, c, i, k̄) + log(α[1])
        p = p - log(-log(rand()))
        p > p̄ && (z[i] = k̄)

        # Update the sufficient statistics
        z̄ != z[i] && suffstats!(s, c, z̄, z[i], i)
    end
end

function update_γ!(
    c::ChainState, 
    s::SuffStats{A, B, C}, 
    pγ0::OffsetArray{Float64, 1, Vector{Float64}}
) where {A, B, C}
    @unpack n, J = s
    @unpack γ, O = c

    # Resample γ[g], given the other γ's
    for g = 2:J
        # log-odds (numerator)
        γ[g] = 1
        suffstats!(s, c)
        log_num = log(pγ0[sum(γ)])
        for k ∈ O, j ∈ (1, g)
            log_num += log_ml(s, j, k)
        end

        # log-odds (denominator)
        γ[g] = 0
        suffstats!(s, c)
        log_den = log(pγ0[sum(γ)])
        for k ∈ O, j ∈ (1)
            log_den += log_ml(s, j, k)
        end

        # log-odds and new γ[g]
        log_odds = log_num - log_den
        γ[g] = rand() <= 1 / (1 + exp(-log_odds))
    end
end

function update_α!(c::ChainState, s::SuffStats{A, B, C}) where {A, B, C}
    @unpack m, N = s
    @unpack α, K = c
    @unpack a0, b0 = m

    ϕ = rand(Beta(α[1] + 1.0, N))
    ψ = (a0 + K[1] - 1.0) / (N * (b0 - log(ϕ))); ψ = ψ / (1 + ψ)
    α[1] = rand(Gamma(a0 + K[1] - (rand() > ψ), 1.0 / (b0 - log(ϕ))))
end

function fit(
    m::Model{A, B},
    y::Vector{C},
    x::Vector{Int};
    iter::Int = 4000, 
    warmup::Int = iter ÷ 2
) where {A, B, C}
    # Initialization
    s = SuffStats(m = m, y = y, x = x)
    c = ChainState(N = s.N, J = s.J)
    pγ1 = zeros(Int, 2^(s.J - 1))
    pγ0 = ph0(s.J - 1, 1.0)

    # Updating
    @fastmath for t ∈ 1:iter
        suffstats!(s, c)
        update_z!(c, s)
        update_α!(c, s)
        update_γ!(c, s, pγ0)
        t > warmup || continue 
        pγ1[γcode(c.γ)] += 1        
    end
    return pγ1
end