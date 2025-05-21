source("npFCBPSFunctional.R")
source("FCBPS.R")
library(RSpectra)
distance2.matrix <- function(X) {
  n <- nrow(X)
  dis2M <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:i) {
      dis2M[i, j] <- sum((X[i, ] - X[j, ])^2)
      dis2M[j, i] <- dis2M[i, j]
    }
  }
  return(dis2M)
}




Gauss.kernel <- function(dis2, h) {
  return(exp(-dis2 / (2 * h^2)))
}






transform <- function(X) {
  Xlim <- apply(X, 2, range)
  Xstd <- matrix(nr = nrow(X), nc = ncol(X))
  for (i in (1:ncol(X))) {
    Xstd[, i] <- (X[, i] - Xlim[1, i]) / diff(Xlim[, i])
  }
  return(list(Xstd = Xstd, Xlim = Xlim))
}



getM <- function(K, thresh.ratio, r = NULL) {
  e <- eigen(K)
  if (is.null(r)) {
    thresh <- e$values[1] * thresh.ratio
    eind <- (abs(e$values) >= thresh)
    r <- sum(eind)
    M <- e$vectors[, eind] %*% sqrt(diag(e$values[eind]))
  } else {
    M <- e$vectors[, 1:r] %*% sqrt(diag(e$values[1:r]))
  }
  return(M)
}

Kw_cal <- function(w, MA, MX) {
  temp <- diag(w) %*% t(KhatriRao(t(MA), t(MX))) - kronecker(MA, t(colMeans(MX)))
  return(temp)
}



eval_obj_grad <- function(w, M1, M2, P, tau) {
  n <- length(w)
  Kw <- P %*% diag(w) %*% M1 - M2
  # KWA <- crossprod(Kw)
  # eigen.out <- eigen(KW)
  eigen.out <- svds(Kw, 1, nu = 0, nv = 1)
  z <- as.numeric(eigen.out$v[, 1])
  v <- M1 %*% z
  colmP <- colMeans(P^2)
  penalty <- tau * sum(w^2 * colmP)
  Ptemp <- sweep(P, 2, v, "*")
  grad <- as.vector(t(2 * Kw %*% z) %*% Ptemp)/n + 2 * tau * w * colmP
  return(list(objective = as.numeric(eigen.out$d[1]^2/n+ penalty), gradient = as.vector(grad)))
}


eval_obj <- function(w, M1, M2, P, tau) {
  n <- length(w)
  # Kw <- cbind(diag(w), -1 * diag(n)) %*% M
  Kw <- P %*% diag(w) %*% M1 - M2
  # KWA <- crossprod(Kw)
  # eigen.out <- eigen(KW)
  eigen.out <- svds(Kw, 1, nu = 0, nv = 1)
  # penalty <- tau * sum((w - 1)^2) / (n^2)
  penalty <- tau * sum(w^2 * colMeans(P^2))
  # obj = eigen.out$values[1] + tau*abs(log(w))
  # obj = eigen.out$values[1] + tau*sum((w-1)^2)
  obj <- as.numeric(eigen.out$d[1]^2)/n + penalty
  return(obj)
}

eval_grad <- function(w, M1, M2, P, tau) {
  n <- length(w)
  Kw <- P %*% diag(w) %*% M1 - M2
  # KWA <- crossprod(Kw)
  # eigen.out <- eigen(KW)
  eigen.out <- svds(Kw, 1, nu = 0, nv = 1)
  z <- as.numeric(eigen.out$v[, 1])
  # v <- Kw %*% z
  v <- M1 %*% z
  Ptemp <- sweep(P, 2, v, "*")
  grad <- as.vector(t(2 * Kw %*% z) %*% Ptemp) /n + 2 * tau * w * colMeans(P^2)
  return(as.vector(grad))
}


barGX_cal <- function(KX, n) {
  barGX <- matrix(0, nrow = n, ncol = n)
  # barGX <- matrix(colMeans(KX), nrow = n, ncol = n, byrow = F)
  rowsum <- rowSums(KX)
  for (i in 1:n) {
    for (j in 1:i) {
      barGX[i, j] <- (rowsum[i] - KX[i, j]) / (n - 1)
      barGX[j, i] <- barGX
    }
  }
  return(barGX)
}


tildeGX_cal <- function(KX, n) {
  tildeGX <- matrix(0, nrow = n, ncol = n)
  sums <- sum(KX)
  for (i in 1:n) {
    for (j in 1:i) {
      tildeGX[i, j] <- sums - sum(KX[i, ]) - sum(KX[j, ]) + KX[i, j]
      tildeGX[j, i] <- tildeGX[i, j]
    }
  }
  tildeGX <- tildeGX / ((n - 1)^2)
  return(tildeGX)
}


cfb_KRR_lam_seq <- function(KA, KAa, Ka, KX, lam_seq, tau = 0, w0 = NULL, lower = NULL, upper = NULL, thresh.ratio = 1e-06, maxit = 2000, xtol_rel = 1e-05, algorithm = "nloptr") {
  n <- dim(KA)[1]
  L <- dim(KAa)[2]

  if (is.null(lower)) {
    lower <- rep(0, n)
  }
  if (is.null(upper)) {
    upper <- rep(Inf, n)
  }
  # n = nrow(Astd)
  if (is.null(w0)) {
    w0 <- rep(1, n)
  }
  # tildeGX <- tildeGX_cal(KX, n)
  tildeGX <- mean(KX)
  # barGX <- barGX_cal(KX, n)
  barGX <- matrix(colMeans(KX), nrow = n, ncol = L, byrow = F)
  GF <- rbind(cbind(KA * KX, KAa * barGX), cbind(t(KAa * (barGX)), Ka * tildeGX))
  # KA = getK_L_prod(Astd)
  # KX = getK_sob_prod(Xstd)

  M <- getM(GF, thresh.ratio)
  M1 <- M[1:n, ]
  M2 <- M[(n + 1):(n + L), ]

  wlist <- matrix(0, nrow = length(lam_seq), ncol = n)
  objlist <- rep(0, length(lam_seq))
  ballist <- rep(0, length(lam_seq))
  for (i in 1:length(lam_seq)) {
    P <- t(solve(KA + lam_seq[i] * diag(n), KAa))
    if (algorithm == "nloptr") {
      res <- nloptr::nloptr(
        x0 = w0, eval_f = eval_obj_grad, lb = lower, ub = upper,
        M1 = M1, M2 = M2, P = P, tau = tau, opts = list(algorithm = "NLOPT_LD_LBFGS", xtol_rel = xtol_rel, maxeval = maxit, check_derivatives = F)
      )
      weights <- res$solution
    }
    if (algorithm == "lbfgsb3c") {
      res <- lbfgsb3c::lbfgsb3c(
        par = w0, fn = eval_obj, gr = eval_grad, lower = lower,
        upper = upper, M1 = M1, M2 = M2, P = P, tau = tau, control = list(
          factr = xtol_rel,
          maxit = maxit
        )
      )
      weights <- res$par
    }
    wlist[i, ] <- weights
    objlist[i] <- res$obj
    w0 <- weights

    Kw <- P %*% diag(weights) %*% M1 - M2
    # KWA <- crossprod(Kw)
    # eigen.out <- eigen(KW)
    eigen.out <- svds(Kw, 1, nu = 0, nv = 1)
    ballist[i] <- eigen.out$d[1]^2/n
  }
  # PPAeval <- gA %*% solve(KA + lambda * diag(n))
  # prodPPAeval <- t(PPAeval) %*% PPAeval



  return(list("wlist" = wlist, "objlist" = objlist, "ballist" = ballist, tau = tau, lam = lam_seq))
}





cfb_KRR <- function(KA, KAa, Ka, KX, lam, tau_seq = exp(seq(log(0.1), log(100), length.out = 10)), w0 = NULL, lower = NULL, upper = NULL, thresh.ratio = 1e-06, maxit = 2000, xtol_rel = 1e-05, algorithm = "nloptr") {
  n <- dim(KA)[1]
  L <- dim(KAa)[2]

  if (is.null(lower)) {
    lower <- rep(0, n)
  }
  if (is.null(upper)) {
    upper <- rep(Inf, n)
  }
  # n = nrow(Astd)
  if (is.null(w0)) {
    w0 <- rep(1, n)
  }
  # tildeGX <- tildeGX_cal(KX, n)
  tildeGX <- mean(KX)
  # barGX <- barGX_cal(KX, n)
  barGX <- matrix(colMeans(KX), nrow = n, ncol = L, byrow = F)
  GF <- rbind(cbind(KA * KX, KAa * barGX), cbind(t(KAa * (barGX)), Ka * tildeGX))
  # KA = getK_L_prod(Astd)
  # KX = getK_sob_prod(Xstd)

  M <- getM(GF, thresh.ratio)
  M1 <- M[1:n, ]
  M2 <- M[(n + 1):(n + L), ]

  wlist <- matrix(0, nrow = length(tau_seq), ncol = n)
  objlist <- rep(0, length(tau_seq))
  ballist <- rep(0, length(tau_seq))
  P <- t(solve(KA + lam * diag(n), KAa))
  for (i in 1:length(tau_seq)) {
    if (algorithm == "nloptr") {
      res <- nloptr::nloptr(
        x0 = w0, eval_f = eval_obj_grad, lb = lower, ub = upper, 
        M1 = M1, M2 = M2, P = P, tau = tau_seq[i], opts = list(algorithm = "NLOPT_LD_LBFGS", xtol_rel = xtol_rel, maxeval = maxit, check_derivatives = F)
      )
      weights <- res$solution
    }
    if (algorithm == "lbfgsb3c") {
      res <- lbfgsb3c::lbfgsb3c(
        par = w0, fn = eval_obj, gr = eval_grad, lower = lower,
        upper = upper, M1 = M1, M2 = M2, P = P, tau = tau_seq[i], control = list(
          factr = xtol_rel,
          maxit = maxit
        )
      )
      weights <- res$par
    }
    wlist[i, ] <- weights
    objlist[i] <- res$obj
    # w0 <- weights

    Kw <- P %*% diag(weights) %*% M1 - M2
    # KWA <- crossprod(Kw)
    # eigen.out <- eigen(KW)
    eigen.out <- svds(Kw, 1, nu = 0, nv = 1)
    ballist[i] <- eigen.out$d[1]^2/n
  }
  # PPAeval <- gA %*% solve(KA + lambda * diag(n))
  # prodPPAeval <- t(PPAeval) %*% PPAeval



  return(list("wlist" = wlist, "objlist" = objlist, "ballist" = ballist, "tau_seq" = tau_seq, "lam" = lam))
}















########## function to output estimated values at evaluation points: use KRR ######
KRR <- function(KA, gA, Y, lam) {
  n <- nrow(KA)
  alpha <- solve(KA + lam * diag(n), Y)
  return(gA %*% alpha)
}


KRR_tune <- function(KA, Y, lam_seq, criterion = "gcv", folds = 5, alpha = 1.4) {
  n <- nrow(KA)
  I <- diag(n)
  valerrors <- rep(0, length(lam_seq))
  if (criterion == "cv") {
    fold_idx <- (1:n) %% folds + 1
  }
  for (i in 1:length(lam_seq)) {
    if (criterion == "gcv" || criterion == "loo") {
      A_mu <- KA %*% solve(KA + lam_seq[i] * I)
    }
    if (criterion == "gcv") {
      # inside <- psych::tr((I - alpha * A_mu) / n)
      inside <- sum(diag((I - alpha * A_mu) / n))
      if (inside < 0) {
        denor <- 0
      } else {
        denor <- (inside)^2
      }
      # denor <- (sum((1 - diag(A_mu)) / n))^2
      numer <- mean(((I - A_mu) %*% Y)^2)
      valerrors[i] <- numer / denor
    }
    if (criterion == "loo") {
      valerrors[i] <- mean(((I - A_mu) %*% Y / (diag(I - A_mu)))^2)
    }
    if (criterion == "cv") {
      for (fold in 1:folds) {
        trainidx <- (fold_idx != fold)
        validx <- (fold_idx == fold)
        Yeval <- KRR(KA[trainidx, trainidx], KA[validx, trainidx], Y[trainidx], lam = lam_seq[i])
        valerrors[i] <- valerrors[i] + mean((Y[validx] - Yeval)^2)
      }
    }
  }
  select_lam <- lam_seq[which.min(valerrors)]
  return(list(select_lam = select_lam, valerrors = valerrors))
}