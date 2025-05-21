# source("perp4_penaltychange.R")
# source("perp4_penalty_center.R")
# library(fdapace)

dis2A <- distance2.matrix(A)
dis2X <- distance2.matrix(X)
##### median heuristic to choose bandwidth in Gaussian kernel
hA <- median(sqrt(dis2A))
hX <- median(sqrt(dis2X))
dis2a <- distance2.matrix(Aeval)

KA <- Gauss.kernel(dis2A, hA)
KX <- Gauss.kernel(dis2X, hX)
Ka <- Gauss.kernel(dis2a, hA)

KAa <- matrix(0, nrow = n, ncol = L)
for (i in 1:neval) {
    for (j in 1:n) {
        dis2 <- sum((Aeval[i, ] - A[j, ])^2)
        KAa[j, i] <- Gauss.kernel(dis2, hA)
    }
}

lam_seq <- exp(seq(log(1e-06), log(1), length.out = 30)) * n


library(mvtnorm)
########## FPCA ##############
# library(fdapace)
# ########## construct Ly and Lt list #########
# L3 <- MakeFPCAInputs(tVec = tgrid, yVec = A)
# FPCA_res <- FPCA(L3$Ly, L3$Lt, list(FVEthreshold = 0.95))
# eigvals <- FPCA_res$lambda
# eigfcs <- FPCA_res$phi #### tgrid times K
# FPCscores <- FPCA_res$xiEst


# L3eval <- MakeFPCAInputs(tVec = tgrid, yVec = Aeval)
# FPCA_res_eval <- predict(FPCA_res, newLy = L3eval$Ly, newLt = L3eval$Lt)
# FPCAscores_eval <- FPCA_res_eval$scores


library(refund)

FPCA_res <- fpca.face(Y= A,  center = TRUE, pve=0.95)
FPCscores <- FPCA_res$scores

FPCA_res <- fpca.face(Y= A, Y.pred = Aeval, center = TRUE, pve=0.95)
FPCAscores_eval <- FPCA_res$scores

source("FCBPS.R")
fcbps <- FCBPS(treat = FPCscores, conf = X)
npfcbps <- npFCBPS.Functional(treat = FPCscores, conf = X, rho = 0.1 / n)

