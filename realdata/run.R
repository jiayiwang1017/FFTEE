
# Load the data
A = readRDS('A.rds')
X = readRDS('X.rds')
Y = readRDS('Y.rds')


source("perp4_penaltychange.R")
finalids <- X$SEQN

###############################
############ get KA and KX
###############################
source("perp4_penaltychange.R")

# dis2A are squared distances between kernel embeddings of A
# saveRDS(dis2A, "./dis2A.rds")
dis2A <- readRDS("./dis2A.rds")
hA <- median(sqrt(dis2A))
KA <- Gauss.kernel(dis2A, hA)
Yidx <- rownames(dis2A) %in% finalids
KA <- KA[Yidx, Yidx]
dim(KA)

Xcon <- transform(as.matrix(X[, c("DMDEDUC2", "RIDAGEYR", "INDFMPIR")]))
X1 <- Xcon$Xstd
dis2X <- distance2.matrix(X1)
hX <- median(sqrt(dis2X))
KXcon <- Gauss.kernel(dis2X, hX)
KX <- KXcon

### construct evaluations on empirical kernel mean embedding functions 

distincA = matrix(c(1:1000), nrow = 1000, ncol = 1)
dis2A <- distance2.matrix(distincA)
h <- median(sqrt(dis2A))
KdistincA <- Gauss.kernel(dis2A, h)
Alist <- list()
for (i in 1:n) {
   Alist[[i]] = as.numeric(A$PAXINTEN[A$SEQN==finalid[i]])
}
agrid <- seq(1, 1000, 5)
Agrid <- matrix(0, nrow = n, ncol = length(agrid))
for (i in 1:n) {
#   print(i)
   a1 = Alist[[i]]
Agrid[i,] = apply(KdistincA[agrid,], 1, function(x, a1) {mean(x[a1])}, a1 = a1)
}
dim(Agrid)
# Agrid <- Agrid[, seq(1, 1000, 5)]
rownames(Agrid) <- finalids





source("perp4_penaltychange.R")

################################################
######## get weights for fcbps and npfcbps#####
################################################
########## FPCA ##############
library(refund)

library(mvtnorm)
########## construct Ly and Lt list #########
FPCA_res <- fpca.face(Y= Agrid,  center = TRUE, pve=0.95)
FPCscores <- FPCA_res$scores
source("FCBPS.R")
fcbps <- FCBPS(treat = FPCscores, conf = X1)
source("npFCBPSFunctional.R")
npfcbps <- npFCBPS.Functional(treat = FPCscores, conf = X1, rho = 0.1 / n)


#################################################
############# get Yhat for other methods ########
#########################
alpha <- 1.4
lam_seq <- exp(seq(log(1e-05), log(0.01), length.out = 10)) * n
lam_seq
out_fc <- function(Yadj, lam_seq, criterion = "loo",  weights = NULL) {
  lam_temp <- KRR_tune(KA, Yadj, lam = lam_seq, criterion = criterion, alpha = alpha)
  lam <- lam_temp$select_lam
  Yhat1 <- KRR(KA, gA = KA, Yadj, lam = lam_temp$select_lam)
  return(list(Yhat = Yhat1, lam = lam))
}
out_nw <- out_fc(Y$BMXBMI, lam_seq)
lam_nw <- out_nw$lam
lam_nw
Yfit_nw <- out_nw$Yhat
summary(Yfit_nw)


# fcbps <- readRDS("./new_subset_results/fcbps_ECDFA_Yidx.rds")
fitdata <- as.data.frame(cbind( Y$BMXBMI, FPCscores))
head(fitdata)

Yfit_fcbps <- out_fc(Yadj = Y$BMXBMI * n * (fcbps$weights), lam_seq, criterion = "loo")
fcbps_lam <- Yfit_fcbps$lam
lm_model <- lm(V1~., data = fitdata, weights = fcbps$weights)
Yfit_fcbps <- predict(lm_model, fitdata)
summary(Yfit_fcbps)


lm_npfcbps <- lm(V1~., data = fitdata, weights = npfcbps$w)
Yfit_npfcbps <- predict(lm_npfcbps, fitdata)
summary(Yfit_npfcbps)



lam_temp <- KRR_tune((KA * KX), ( Y$BMXBMI), lam = lam_seq, criterion = "loo", alpha = alpha)
print(lam_temp)
gA <- KA * matrix(colMeans(KX), nrow = n, ncol = n, byrow = T)
Yfit_reg <- KRR((KA * KX), gA = gA, ( Y$BMXBMI), lam = lam_temp$select_lam)
Yfit <- KRR((KA * KX), gA = (KA * KX), ( Y$BMXBMI), lam = lam_temp$select_lam)
summary(Yfit_reg)


train_data <- as.data.frame(cbind(Y$BMXBMI, FPCscores, Xcon$Xstd))
names(train_data) <- c("Y", "A1", "A2", "X1", "X2", "X3")
summary(train_data)
lm_fit <- lm(Y ~ ., data = train_data)

TX <- mean(Xcon$Xstd %*% lm_fit$coefficients[c("X1", "X2", "X3")] )
coef <- as.numeric(lm_fit$coefficients[!names(lm_fit$coefficients) %in% c("(Intercept)","X1", "X2", "X3")])
Yfit_flm <- lm_fit$coefficients[1] + as.matrix(FPCscores) %*% coef + TX
summary(Yfit_flm)



X1 <- Xcon$Xstd[,1]
X2 <- Xcon$Xstd[,2]
X3 <- Xcon$Xstd[,3]
BMI <- as.numeric(Y$BMXBMI)
gamfit <- pfr(BMI~ af(Agrid, k=c(6,6), m =list(c(2,3), c(2,3)), bs="ps",
presmooth="fpca.face") + s(X1, bs = "cr", k = 5) + s(X2, bs = "cr", k =5) + s(X3, bs = "cr", k=5) )
get_effect <- function(a){
Atemp <- matrix(a, nrow = nrow(Agrid), ncol = ncol(Agrid), byrow = TRUE)
newdata <- list(Agrid = Atemp, X1 = X1, X2 = X2, X3 = X3)
pred <- predict(gamfit, newdata=newdata)
return(mean(pred))
}

get_mult_effect <- function(Aeval){
    Yeffect <- rep(0, nrow(Aeval))
    for (i in 1:nrow(Aeval)){
        a = Aeval[i,]
        Yeffect[i] = get_effect(a)
    }
    return(Yeffect)
}

Yfit_gam <- get_mult_effect(Agrid)
summary(Yfit_gam)


######################################################
################### get Yhat for cfb  #################
######################################################
tau_seq <- c(0, (exp(seq(log(1e-03), log(1), length.out = 10))))
out <- cfb_KRR(KA, KA, KA, KX, lam = fcbps_lam, tau_seq = tau_seq, maxit = 6000, xtol_rel = 1e-8)

comp_errors <- function(Y, Yfit, Yfitadj, lam, tau_idx){
        comp1 =   KRR(KA, gA = KA, Yfit * out$wlist[tau_idx,], lam = lam) - Yfitadj
        comp2 = KRR(KA, gA = KA, (Y-Yfit) * out$wlist[tau_idx,], lam = lam) 
        return(mean(comp1^2) + mean(comp2^2))
    }
test_errors <- rep(0, length(out$tau_seq))
for(i in 1:length(out$tau_seq)){
    test_errors[i] = comp_errors(Y$BMXBMI, Yfit, Yfit_reg, out$lam, tau_idx = i)
}
# test_errors
select_tau_idx <- which.min(test_errors)
summary(out$wlist[select_tau_idx,])
Yfit_cfb <- KRR(KA, gA = KA, Y$BMXBMI * (out$wlist[select_tau_idx,]), lam = out$lam)
summary(Yfit_cfb)


######## Final Fitted Results ##########
Yfit_final <- cbind(Yfit_nw, Yfit_fcbps, Yfit_npfcbps, Yfit_reg, Yfit_flm, Yfit_gam, Yfit_cfb)
Yfit_final <- as.data.frame(Yfit_final)
names(Yfit_final) <- c("nw", "fcbps", "npfcbps", "reg", "flm", "fgam", "cfb")
