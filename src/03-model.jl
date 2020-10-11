"""
    Hyperparameters of our MANOVA-DDP model, following the paper notation.
"""
@with_kw struct Model{A, B}
    D::Int
    r0::Int = 1
    ν0::Int = D + 2
    u0::A = zeros(D)
    S0::B = cholesky(collect(I(D)))
    a0::Float64 = 1.0
    b0::Float64 = 1.0
    ζ0::Float64 = 1.0
end