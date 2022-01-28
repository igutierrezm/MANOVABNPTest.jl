# Fit our model to the data in data/allison1962.csv and plot the results

#==============================================================================#
# Part I - Fitting
#==============================================================================#

# Load the relevant R packages
library(dplyr)
library(ggplot2)
library(readr)
library(JuliaConnectoR) # includes juliaEval, juliaImport, juliaGet

# Install the relevant Julia packages
juliaEval('Pkg.add(url = "https://github.com/igutierrezm/MANOVABNPTest.jl")')

# Load the relevant Julia packages from our R session
MANOVABNPTest <- juliaImport('MANOVABNPTest')
Random <- juliaImport('Random')

# Load the data
df <- readr::read_csv("data/allison1962.csv")

# Compute y, x and grid
x <- as.integer(df[[1]])
y <- as.matrix(df[, 2:3])
y[, 1] <- scale(y[, 1])
y[, 2] <- scale(y[, 2])
grid <- seq(-4, 4, length.out = 100)

# Fix the seed
rng <- Random$MersenneTwister(1L)

# Fit the model
out <- MANOVABNPTest$train(y, x, grid, iter = 10000L, rng = rng)
 
# Extract the posterior of the key objects
posterior_g <- data.frame(out$hypotheses) # p(gamma \mid y, x)
posterior_f <- data.frame(out$densities)  # E[f \mid y, x]

#==============================================================================#
# Part II - Plotting
#==============================================================================#

# Set the group labels
lbls <- c(
    "control", 
    "treatment-1", 
    "treatment-2",
    "treatment-3"
)

# Compute the plot points
pts <- 
    readr::read_csv("data/allison1962.csv") %>%
    dplyr::mutate(
        group = factor(group, labels = lbls),
        y1 = (y1 - mean(y1)) / sd(y1),
        y2 = (y2 - mean(y2)) / sd(y2)
    )

# Compute the plot contours
cntrs <-
    posterior_f %>%
    dplyr::mutate(group = factor(j, labels = lbls))

# Create the final plot
p <- 
    cntrs %>%
    ggplot(aes(colour = group)) +
    geom_contour(aes(x = y1, y = y2, z = f), binwidth = 0.05, size = 0.2) +
    geom_point(data = pts, aes(y1, y2, colour = group)) +
    labs(
        x = "number of bacilli inhaled per tubercle formed (z-score)",
        y = "tubercle size (z-score)",
        colour = "group"
    ) +
    theme_linedraw() +
    theme(
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position  = "top"
    )
p
