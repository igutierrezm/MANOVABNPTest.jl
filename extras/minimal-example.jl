# Activate the environment
using Pkg
Pkg.activate("extras")

# Install MANOVABNPTest and Random
Pkg.add(url = "https://github.com/igutierrezm/MANOVABNPTest.jl")
Pkg.add("Random")

# Load the Julia packages
using MANOVABNPTest
using Random

# Simulate (y, x, grid)
begin 
    rng = MersenneTwister(1)
    grid = collect(-3.0:0.1:3.0)
    x = 1 .+ collect(0:99) .% 4
    y = randn(rng, 100, 2)
end;

# Fit the model to our simulated data
out = train(y, x, grid; rng);

# Print the posterior hypothesis distribution
first(out.hypotheses, 3)
# 5×2 DataFrame
#  Row │ hypothesis  prob    
#      │ String      Float64 
# ─────┼─────────────────────
#    1 │ [0, 0, 0]    0.9835
#    2 │ [1, 0, 0]    0.015
#    3 │ [0, 1, 0]    0.0005

# Print the posterior predictive density
first(out.densities, 3)
# 3×6 DataFrame
#  Row │ j      var1   var2   y1       y2       f           
#      │ Int64  Int64  Int64  Float64  Float64  Float64     
# ─────┼────────────────────────────────────────────────────
#    1 │     1      1      2     -3.0     -3.0  0.000136175
#    2 │     1      1      2     -2.9     -3.0  0.000161726
#    3 │     1      1      2     -2.8     -3.0  0.000191867
