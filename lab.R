library(tidyverse)
library(Matrix)
library(spBayes)


n <- 400
coords <- expand.grid(x <- seq(0, 1, length.out = sqrt(n)), x)
plot(coords, pch = 19, cex = 0.8)

## the model has difficulty separating phi and sigma_sq as resolution / number
## of points increases, so their samples are inversely correlated
sigmasq <- 2
phi <- 0.5

C <- sigmasq * exp(-phi * as.matrix(dist(coords)))

# expoonential covariance for a GP = Ornstein-Uhlenbeck process in 2d
# on a regular grid of location --> block Toeplitz with Toeplitz blocks matrix
C %>% as("dgCMatrix") %>% image()

# Simulate the latent process
# chol gives you upper triangular
L <- t(chol(C))
L %>% as("dgCMatrix") %>% image()

w <- L %*% rnorm(n, 0, 1)
df <- data.frame(coords, w)

ggplot(df, aes(x = Var1, Var2, fill = w)) + geom_raster() +
    scale_fill_viridis_c()


# Issue here - the model cannot identify a difference between w and beta_0
p <- 3
X <- cbind(1, matrix(rnorm(n * (p - 1)), nrow = n))
beta <- c(2.5, -1.5, 4.0)


eps <- 0.1^(0.5) * rnorm(n)

y <- X %*% beta + w + eps


df <- data.frame(coords, w, y)
ggplot(df, aes(Var1, Var2, fill = y)) + geom_raster() +
    scale_fill_viridis_c()


n.samples <- 2000
starting <- list("phi" = 3 / 0.5,
                 "sigma.sq" = 50,
                 tau.sq = 1)
tuning <- list("phi" = 0.1,
               "sigma.sq" = 0.1,
               "tau.sq" = 0.1)
priors.1 <- list(
    "beta.Norm" = list(rep(0, p), diag(1000, p)),
    "phi.Unif" = c(0.1, 3 / 0.1),
    "sigma.sq.IG" = c(2, 2),
    "tau.sq.IG" = c(2, 0.1)
)

cov.model <- "exponential"
n.report <- 5000
n.samples <- 2000
verbose = TRUE


# discussion on centering constraints within MCMC:
# bayesian backfitting
m.1 <- spLM(
    y ~ X - 1,
    coords = as.matrix(coords),
    starting = starting,
    tuning = tuning,
    prior = priors.1,
    cov.model = cov.model,
    n.samples = n.samples,
    verbose=verbose,
    n.report=n.report
)

burn.in <- 0.5*n.samples
m.1 <- spRecover(m.1, start=burn.in, verbose=FALSE)


colMeans(m.1$p.beta.recover.samples)

plot(m.1$p.beta.recover.samples[, 1], type = "l")

# Notice that there are issues with identifiability in several parameters,
# most notably around the intercept and sigma^2/phi parameters
