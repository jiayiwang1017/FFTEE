phi1 <- function(x, k) {
    return(sqrt(2) * sin(2 * pi * k * x))
}

phi2 <- function(x, k) {
    return(sqrt(2) * cos(2 * pi * k * x))
}

# phi <- function(x, k) {
# if(k %%2 == 1){return(phi1(x,k))}else {
#   return(phi2(x,k))
# }
# }


phi <- function(x, k) {
    if (k %% 2 == 1) {
        return(phi1(x, (k + 1) / 2))
    } else {
        return(phi2(x, k / 2))
    }
}


grid_size <- 100
tgrid <- seq(from = 0, to = 1, length.out = grid_size)

Z1 <- rnorm(n)
Z2 <- rnorm(n)
Z3 <- rnorm(n)
Z4 <- rnorm(n)
Z <- cbind(Z1, Z2, Z3, Z4)


Z1_ <- rnorm(n)
Z2_ <- rnorm(n)
Z3_ <- rnorm(n)
Z4_ <- rnorm(n)
Z_ <- cbind(Z1_, Z2_, Z3_, Z4_)

##### FPC
A1 <- 4 * Z1 + rnorm(n)
A2 <- 2 * sqrt(3) * Z2 + rnorm(n)
A3 <- 2 * sqrt(2) * Z3 + rnorm(n)
A4 <- 2 * Z4 + rnorm(n)
AF <- cbind(A1, A2, A3, A4)
K <- 4
trueeigvals <- c(17, 13, 9, 5)
A <- matrix(0, nrow = n, ncol = grid_size)

for (k in 1:K) {
    A <- A + AF[, k] %o% (phi(tgrid, k))
}


##### Covariates (confounders)
X <- cbind(Z1, Z2, Z3, Z4)

#### true causal effect
mu <- 2 * sqrt(2) * sin(2 * pi * tgrid) + sqrt(2) * cos(2 * pi * tgrid) + sqrt(2) * sin(4 * pi * tgrid) / 2 + sqrt(2) * cos(4 * pi * tgrid) / 2


myfun <- function(x, mu) {
    return(mean(x * mu))
}

Yeffect <- apply(A, 1, myfun, mu = mu)
######### random evaluation points ########
neval <- 100
L <- neval
library(MASS, help, pos = 2, lib.loc = NULL)
eval_points <- mvrnorm(n = neval, mu = rep(0, length(trueeigvals)), Sigma = diag(trueeigvals))

Aeval <- matrix(0, nrow = neval, ncol = grid_size)
for (k in 1:K) {
    Aeval <- Aeval + eval_points[, k] %o% (phi(tgrid, k))
}

Yeval <- apply(Aeval, 1, myfun, mu = mu)