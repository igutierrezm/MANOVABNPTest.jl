# Activate a Julia environment
library(JuliaConnectoR)
Pkg <- juliaImport('Pkg')
Pkg$activate("extras")

# Install MANOVABNPTest and Random
Pkg$add(url = "https://github.com/igutierrezm/MANOVABNPTest.jl")
Pkg$add("Random")

# Load MANOVABNPTest and Random
MANOVABNPTest <- juliaImport('MANOVABNPTest')
Random <- juliaImport('Random')

# Simulate (y, x, grid)
set.seed(1)
y <- matrix(rnorm(200), 100, 2)
x <- as.integer(1 + 0:99 %% 4)
grid <- seq(-3, 3, by = 0.1)

# Fit the model to our simulated data
rng <- Random$MersenneTwister(1L)
out <- MANOVABNPTest$train(y, x, grid, rng = rng)

# Collect the results
posterior_g <- data.frame(out$hypotheses) # p(gamma | y, x)
posterior_f <- data.frame(out$densities)  # p(y* | y, x)

# Print the posterior hypothesis distribution
head(posterior_g, 3)
#   hypothesis   prob
# 1  [0, 0, 0] 0.7715
# 2  [1, 0, 0] 0.0055
# 3  [0, 1, 0] 0.2120

# Print the posterior predictive density
head(posterior_f, 3)
#   j var1 var2   y1 y2            f
# 1 1    1    2 -3.0 -3 2.802552e-05
# 2 1    1    2 -2.9 -3 3.213636e-05
# 3 1    1    2 -2.8 -3 3.708808e-05
