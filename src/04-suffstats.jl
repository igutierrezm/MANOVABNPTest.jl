"""
    Sufficient statistics of our MANOVA-DDP model, 
    see Bouchard-Côté et al. (2017, eqs. 28-29)
"""
@with_kw struct SuffStats{A, B, C}
    m::Model{A, B}
    y::Vector{C}
    x::Vector{Int}
    K̄::Int                     = 64
    N::Int                     = length(y)
    J::Int                     = maximum(x)
    D::Int                     = length(y[1])
    n::Vector{Int}             = zeros(Int, K̄)        
    r::Matrix{Int}             = [m.r0 for j = 1:J, k ∈ 1:K̄]        
    ν::Matrix{Int}             = [m.ν0 for j = 1:J, k ∈ 1:K̄]        
    u::Matrix{A}               = [deepcopy(m.u0) for j = 1:J, k ∈ 1:K̄]    
    S::Matrix{B}               = [deepcopy(m.S0) for j = 1:J, k ∈ 1:K̄]
    Σ::Matrix{Matrix{Float64}} = [zeros(D, D) for j = 1:J, k ∈ 1:K̄]
    μ::Matrix{Vector{Float64}} = [zeros(D) for j = 1:J, k ∈ 1:K̄]
end

"""
    Recompute the sufficient statistics from scratch
"""
function suffstats0!(s::SuffStats{A, B, C}, c::ChainState) where {A, B, C}
    @unpack m, y, x, K̄, N, J, n, r, ν, u, S = s
    @unpack r0, ν0, u0, S0 = m
    @unpack z, γ, O = c
    idx = Bool.(ones(N))
    @inbounds for j = 1:J, k ∈ 1:K̄
        n[k] = 0
        @. idx = (x^γ[x] == j) & (z == k)
        if any(idx)
            ysub = y[idx]
            njk = length(ysub)
            ym = mean(ysub)
            yv = cov(ysub, corrected = false)
            ν[j, k] = ν0 + njk
            r[j, k] = r0 + njk
            u[j, k] = (r0 * u0 + njk * ym) / r[j, k]
            S[j, k] = 
                (Matrix(S0) + reinterpret(Float64, njk * yv  + njk * r0 * (ym - u0) * (ym - u0)' / r[j, k])) |>
                x -> Symmetric(x) |>
                x -> cholesky(x)
        else
            ν[j, k] = ν0
            r[j, k] = r0
            u[j, k] .= u0
            S[j, k].factors .= S0.factors
        end
    end
    @inbounds @fastmath for i ∈ 1:N 
        n[z[i]] += 1 
    end
    return s
end

"""
    Recompute the sufficient statistics from scratch
"""
function suffstats!(s::SuffStats{A, B, C}, c::ChainState) where {A, B, C}
    @unpack m, y, x, K̄, N, J, n, r, ν, u, S = s
    @unpack r0, ν0, u0, S0 = m
    @unpack z, γ = c

    @inbounds @fastmath for j ∈ 1:J, k ∈ 1:K̄
        n[k] = 0
        ν[j, k] = ν0
        r[j, k] = r0
        u[j, k] .= u0
        S[j, k].factors .= S0.factors
    end

    @inbounds @fastmath for i ∈ 1:N
        j = x[i]^γ[x[i]]
        k = z[i]
        n[k] += 1
        νi = ν[j, k] += 1
        ri = r[j, k] += 1
        u[j, k] .= ((ri - 1) * u[j, k] + y[i]) / ri
        u[j, K̄] .= (y[i] - u[j, k]) * √(ri / (ri - 1))
        lowrankupdate!(S[j, k], u[j, K̄])
    end
    return s
end

"""
    Recompute the sufficient statistics when observation `i0` moves from 
    cluster `k0` to `k1`, see Bouchard-Côté et al. (2017, eq. 30)
"""
function suffstats!(
    s::SuffStats{A, B, C},
    c::ChainState,
    k0::Int,
    k1::Int,
    i0::Int
) where {A, B, C}
    @unpack m, y, x, n, r, ν, u, S, K̄ = s
    @unpack r0, ν0, u0, S0 = m
    @unpack γ, z, U, O, K = c
    xi = x[i0]^γ[x[i0]]
    yi = y[i0]

    ##### Modify group/cluster j/k1 ############################################

    # Open the cluster, if pertinent
    if n[k1] == 0 
        push!(O, k1)
        pop!(U, k1)
        K[1] += 1
    end

    # Update the sufficient stats
    n[k1] += 1
    uk = u[xi, k1]
    rk = r[xi, k1] += 1
    νk = ν[xi, k1] += 1
    uk .= ((rk - 1) * uk + yi) / rk
    u[xi, K̄] .= (yi - uk) * √(rk / (rk - 1))
    lowrankupdate!(S[xi, k1], u[xi, K̄])

    ##### Modify group/cluster j/k0 ############################################

    # Update the sufficient stats
    n[k0] -= 1
    uk = u[xi, k0]
    rk = r[xi, k0] -= 1
    νk = ν[xi, k0] -= 1
    u[xi, K̄] .= (yi - uk) * √((rk + 1) / rk)
    lowrankdowndate!(S[xi, k0], u[xi, K̄])
    uk .= ((rk + 1) * uk - yi) / rk

    # Close the cluster, if pertinent
    if n[k0] == 0 
        push!(U, k0)
        pop!(O, k0)
        K[1] -= 1
    end
    return s
end
