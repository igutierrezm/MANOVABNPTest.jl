# MANOVABNPTest.jl

A Julia package that implements the Bayesian nonparametric (BNP) MANOVA model described in 
Gutiérrez et al. (2022).

## Installation

Install with the Julia package manager Pkg:

```julia
# Press ']' to enter the Pkg REPL mode.
pkg> add https://github.com/igutierrezm/MANOVABNPTest.jl
```
or
```julia
julia> using Pkg; 
julia> Pkg.add(url = "https://github.com/igutierrezm/MANOVABNPTest.jl")
```
We recommend you create a Pkg environment for each project you use JuMP for, instead of adding lots of packages to the global environment. The [Pkg manager documentation](https://julialang.github.io/Pkg.jl/v1/environments/) has more information on this topic. For example, if you want to install `MANOVABNPTest.jl` to use it in a project located at `<dir>`, you should install `MANOVABNPTest.jl` as follows:
```julia
julia> using Pkg; 
julia> Pkg.activate("<dir>")
julia> Pkg.add(url = "https://github.com/igutierrezm/MANOVABNPTest.jl")
```

## Getting started

After installation, you can execute the hypothesis test by calling `train()`. This function has 3 main arguments: `y` (an `N x D` matrix of outcomes), `x` (an `N`-dimensional vector of group labels) and `grid` (an `M`-dimensional vector indicating `M` grid points for plotting purposes). You can also add a random number generator `rng` if desired. Here is a minimal example:
```julia
# Load the relevant datasets
using MANOVABNPTest
using Random

# Simulate (y, x, grid)
rng = MersenneTwister(1)
grid = collect(-3.0:0.1:3.0)
x = 1 .+ collect(0:99) .% 4
y = randn(rng, 100, 2)

# Fit the model
out = train(y, x, grid; rng);
```
The result is a named tuple with 2 `DataFrame`s: `hypotheses` (containing the posterior probability of each hypothesis) and `densities` (containing the posterior predictive distribution over a cartesian product of the selected grid of points). Both dataframes are already in tidy format, so they are easy to use in combination with any ploting library that explotes the grammar of graphics, such as [`AlgebraOfGraphics.jl`](http://juliaplots.org/AlgebraOfGraphics.jl/dev/).

## Using MANOVABNP.jl from R

R users can also use this package thanks to the R package `JuliaConnectoR`. Verify, install Julia.

> **Tip:** We recommend installing Julia following these [platform-specific instructions](https://julialang.org/downloads/platform/).

Next, make sure that JuliaConnectoR can call Julia:

> **Tip:** The easiest way to ensure this is by adding Julia's `bin/` folder to your `PATH` environment variable. This is also explained in detail in Julia's [platform-specific instructions](https://julialang.org/downloads/platform/). It is important to understand that installing Julia is not enough, it must also be reachable from R.

Next, install `JuliaConnectoR` as usual:
```R
install.packages("JuliaConnectoR")
```
Next, create/activate a Julia environment. For example, if your project is located as `<dir>`, you should create the environment as follows:
```R
library(JuliaConnectoR)
Pkg <- juliaImport('Pkg')
Pkg$activate("<dir>")
```
Next, install `MANOVABNPTest.jl` as follows:
```R
Pkg$add(url = "https://github.com/igutierrezm/MANOVABNPTest.jl")
```
Next, import the Julia packages `MANOVABNPTest` and `Random` as follows:
```R
MANOVABNPTest <- juliaImport('MANOVABNPTest')
Random <- juliaImport('Random')
```
After this setup, we can execute our hypothesis test in R just as in Julia. Here is a minimal example:
```R
set.seed(1)
y <- matrix(rnorm(200), 100, 2)
x <- as.integer(1 + 0:99 %% 4)
grid <- seq(-3, 3, by = 0.1)
rng <- Random$MersenneTwister(1L)
out <- MANOVABNPTest$train(y, x, grid, rng = rng)
```
Note that the result is a pointer to a Julia object. To collect the results, we can proceed as follows:
```R
posterior_g <- data.frame(out$hypotheses) # p(gamma | y, x)
posterior_f <- data.frame(out$densities)  # p(y* | y, x)
```

Here is a glimpse the collected results
```R
head(posterior_g, 3)
#   hypothesis   prob
# 1  [0, 0, 0] 0.7715
# 2  [1, 0, 0] 0.0055
# 3  [0, 1, 0] 0.2120

head(posterior_f, 3)
#   j var1 var2   y1 y2            f
# 1 1    1    2 -3.0 -3 2.802552e-05
# 2 1    1    2 -2.9 -3 3.213636e-05
# 3 1    1    2 -2.8 -3 3.708808e-05
```

As you can see, the results are already in tidy format, so they are easy to use in combination with any ploting library that explotes the grammar of graphics, such as [`ggplot2`](https://ggplot2.tidyverse.org/). However, the columns of inside the dataframe `posterior_f` are not easy to interpret. Here is the definition of each variable:

- `f`: The value of posterior predictive density when y[`var1`] == `y1`, y[`var2`] == `y2` and x == `j`.
- `j`: The group indicator.
- `var1`, `var2`: Which outcomes are considered in this row.
- `y1`: The value of the variable with ID `var1` considered in this row.
- `y2`: The value of the variable with ID `var2` considered in this row.

Using this information, you can then plot each posible bivariate density.

## Gallery

More elaborate examples are stored in `extras/`. 

## Notes

1. Currently, this package is a prototype. I'm planning to convert this prototype into a fully fledged package (including a proper documentation) during the following 6 months.

2. Currently, the package has been tested under Julia 1.5.3 and R 4.0.0. For full reproducibility, consider running this example in a Docker container generated from the Docker image generated by the file `Dockerfile`.
