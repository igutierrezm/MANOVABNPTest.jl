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
We recommed to activate a local environment before installing this package. For example, if you want to install `MANOVABNPTest.jl` to use it in a project located at `<dir>`, you should install `MANOVABNPTest.jl` as follows:
```julia
julia> using Pkg; 
julia> Pkg.activate("<dir>")
julia> Pkg.add("MANOVABNPTest")
```

## Getting started

After installation, you can execute the hypothesis test by calling `train()`. This function has 3 main arguments: `y` (an `N x D` matrix), `x` (an `N`-dimensional vector of group labels) and `grid` (an `M`-dimensional grid indicating `M` grid points for the purposes of plotting). You can also add a random number generator `rng` if desired. Here is a minimal example:
```
# Load the relevant datasets
using MANOVABNPTest
using Random

# Simulate (y, x, grid)
y = randn(100, 2);
x = 1 .+ collect(0:99) .% 4;
grid = collect(-3.0:0.1:3.0);

# Fit the model
rng = MersenneTwister(1);
out = train(y, x, grid; rng);
```

The result is a named tuple with 2 `DataFrame`s: `hypotheses` (containing the posterior probability of each hypothesis) and `densities` (containing the posterior predictive distribution over a cartesian product of the selected grid of points). Both dataframes are already in tidy format, so they are easy to use in combination with any ploting library that explotes the grammar of graphics, such as `AlgebraOfGraphics.jl`.

## Using MANOVABNP.jl from R

R users can also epxloit this package using the R package `JuliaConnectoR`. First, install R and Julia. Next, install `JuliaConnectoR` as follows:
```R
R> install.packages("JuliaConnectoR")
```
Next, activate a Julia environment. For example, if your project is located as `<dir>`, you should create the environment as follows:
```R
# Activate a Julia environment
R> library(JuliaConnectoR)
R> Pkg <- juliaImport('Pkg')
R> Pkg$activate("<dir>")
```
Next, install `MANOVABNPTest.jl` as follows:
```R
R> juliaEval('Pkg.add(url = "https://github.com/igutierrezm/MANOVABNPTest.jl")')
```
Next, import the Julia packages `MANOVABNPTest` and `Random` as follows:
```R
R> MANOVABNPTest <- juliaImport('MANOVABNPTest')
R> Random <- juliaImport('Random')
```
After this setup, we can execute our hypothesis test in R just as in Julia. Here is a minimal example:
```R
y <- matrix(rnorm(200), 100, 2)
x <- as.integer(1 + 0:99 %% 4)
grid <- seq(-3, 3, by = 0.1)
rng <- Random$MersenneTwister(1L)
out <- MANOVABNPTest$train(y, x, grid, rng = rng)
```
Note that the result is a pointer to a Julia object. To collect the results, we can proceed as follows:
```R
R> posterior_g <- data.frame(out$hypotheses) # p(gamma | y, x)
R> posterior_f <- data.frame(out$densities)  # p(y* | y, x)
```

## Gallery

More elaborate examples are stored in `extras/`. 

## Notes

1. Currently, this package is a prototype. I'm plannig to convert this prototype into a fully fledged package (including a proper documentation) during the following 6 months.

2. Currently, the package has been tested under Julia 1.5.3 and R 4.0.4. For full reproducibility, consider running this example in a Docker container generated from this [Docker image](https://hub.docker.com/r/igutierrez1988/manovabnptest-example).
