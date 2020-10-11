"""
    Compute ln p(y[i] | y[-i], z[-i], z[i] = k), given `i` and `k`,
    see Bouchard-Côté et al. (2017, eq. 32) and Sudderth (2006, eq. 2.163)
"""
function log_pl(
    s::SuffStats{A, B, C}, 
    c::ChainState, 
    i::Int, 
    k::Int
) where {A, B, C}
    @unpack y, x, K̄, D, ν, r, u, S = s
    @unpack γ, z = c
    yi = y[i]
    zi = z[i]
    j  = x[i]^γ[x[i]]

    if zi == k
        r1 = r[j, k];   r0 = r1 - 1
        ν1 = ν[j, k];   ν0 = ν1 - 1    
        S1 = S[j, k]
        S0 = S[j, K̄]
        S0.factors .= S1.factors
        u[j, K̄] .= (yi - u[j, k]) * √(r1 / r0)
        lowrankdowndate!(S0, u[j, K̄])
    else
        r0 = r[j, k];   r1 = r0 + 1
        ν0 = ν[j, k];   ν1 = ν0 + 1
        S0 = S[j, k]
        S1 = S[j, K̄]
        S1.factors .= S0.factors
        u[j, K̄] .= (r0 * u[j, k] + yi) / r1
        u[j, K̄] .= (yi - u[j, K̄]) * √(r1 / r0)
        lowrankupdate!(S1, u[j, K̄])
    end
    
    val =
        0.5ν0 * logdet(S0) -
        0.5ν1 * logdet(S1) +
        logmvgamma(D, ν1 / 2) -
        logmvgamma(D, ν0 / 2) +
        0.5D * log(r0 / r1) -
        0.5D * log(π)
end

"""
    Compute ln p({y[i]: x[i] = j, z[i] = k} | x, z), given `j` and `k`,
    see Bouchard-Côté et al. (2017, eq. 31)
"""
function log_ml(s::SuffStats{A, B, C}, j::Int, k::Int) where {A, B, C}
    @unpack m, D, r, ν, S = s
    @unpack r0, ν0, S0 = m
    r1 = r[j, k]
    ν1 = ν[j, k]
    S1 = S[j, k]
    ll =
        0.5ν0 * logdet(S0) -
        0.5ν1 * logdet(S1) +
        logmvgamma(D, ν1 / 2) -
        logmvgamma(D, ν0 / 2) +
        0.5D * log(r0 / r1) -
        0.5D * log(π) * (r1 - r0)
end
