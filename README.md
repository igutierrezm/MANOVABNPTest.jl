# MANOVABNPTest.jl

A Julia package that implements the (BNP) MANOVA model described in 
Gutiérrez et al. (2021). See the documentation for details and the getting started vignette for a working example.

## Installation

Install with the Julia package manager Pkg:

```julia
# Press ']' to enter the Pkg REPL mode.
pkg> add PolyaGammaSamplers
```
or
```julia
julia> using Pkg; 
julia> Pkg.add("MANOVABNPTest")
```

## Getting started

After installation, you can execute the hypothesis test by calling `train()`. This function has 3 main arguments: `y` (an $N \times D$ matrix), $x$ (an $N$-dimensional vector of group labels) and `grid` (an $M$-dimensional grid indicating $M$ grid points for the purposes of plotting). Here is a minimal example:
```
using MANOVABNPTest, Random
N = 200
y = randn(200, 2)
x = 1 .+ collect(0:199) .% 4
rng = MersenneTwister(1)
grid = collect(-4.0:0.1:4.0)
out = train(y, x, grid; rng)
```

The result is a named tuple with 2 `DataFrame`s: `hypotheses` (containing the posterior probability of each hypothesis) and `densities` (containing the posterior predictive distribution over a cartesian product of the selected grid of points). Both dataframes are already in tidy format, so they are easy to use in combination with any ploting library that explotes the grammar of graphics, such as `AlgebraOfGraphics.jl`.

## Using MANOVABNP.jl from R

R users can also epxloit this package using the R package `JuliaConnectoR`. First, install R and Julia. Next, install `JuliaConnectoR` as follows:
```R
R> install.packages("JuliaConnectoR")
```
Next, install `MANOVABNPTest.jl` as follows:
```R
R> library(JuliaConnectoR)
R> juliaEval('Pkg.add(url = "https://github.com/igutierrezm/MANOVABNPTest.jl")')
```
Next, import the Julia packages `MANOVABNPTest` and `Random` (neccesary to fix the seed) as follows:
```R
R> MANOVABNPTest <- juliaImport('MANOVABNPTest')
R> Random <- juliaImport('Random')
```
After this setup, we can execute our hypothesis test in R just as in Julia. Here is a minimal example:
```R
y <- matrix(rnorm(400), 200, 2)
x <- as.integer(1 + 0:199 %% 4)
grid <- seq(-4, 4, by = 0.1)
rng <- Random$MersenneTwister(1L)
out <- MANOVABNPTest$train(y, x, grid, rng = rng)
```
Note that the result is a pointer to a Julia object. To collect the results, we can proceed as follows:
```R
R> posterior_g <- data.frame(out$hypotheses) # p(gamma \mid y, x)
R> posterior_f <- data.frame(out$densities)  # E[f \mid y, x]
```

## Notes

1. Currently, this package is a prototype. I'm plannig to convert this prototype into a fully fledged package (including a proper documentation) during the following 6 months.

2. Currently, the package has been tested under Julia 1.5.3 and R 4.0.4. For full reproducibility, consider running this example in a Docker container generated from the image \url{https://hub.docker.com/r/igutierrez1988/manovabnptest-example}.
