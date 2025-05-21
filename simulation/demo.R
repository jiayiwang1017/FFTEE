### sample size n
### functional treatment A
### confounders X
### noise std sigma

### causal effect Yeval(on evaluation points) 



n <- 200
sigma <- 1

### generate data
source("data_generate_linear.R")
X1 <- X[,1]
X2 <- X[,2]
X3 <- X[,3]
X4 <- X[,4]
### outcome model 
YX <-  (Z2 * Z1^2 + Z4^2 * sin(2 * Z3))
Y <- 15 * YX + Yeffect + rnorm(n, 0, sd = sigma)



### compute weights for fcbps and npfcbps
source("perp4_penaltychange.R")
source("cal_weights.R")


#### construct gram matrices
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

#### get FPCA
library(refund)

FPCA_res <- fpca.face(Y= A,  center = TRUE, pve=0.95)
FPCscores <- FPCA_res$scores

FPCA_res <- fpca.face(Y= A, Y.pred = Aeval, center = TRUE, pve=0.95)
FPCAscores_eval <- FPCA_res$scores


#### nw estimator 
lam_seq <- exp(seq(log(1e-06), log(1), length.out = 30)) * n
lam_nw <- KRR_tune(KA, Y, lam_seq, criterion = "loo", alpha = 1)
Yevalhat <- KRR(KA, gA = t(KAa), Y, lam = lam_nw$select_lam)
nw_mse <- mean((Yeval - Yevalhat)^2)
print(nw_mse)



### fcbps estimator
FPCscores <- as.data.frame(FPCscores)
train_data <- as.data.frame(cbind(Y, FPCscores))
summary(train_data)
fcbps_lm <- lm(Y ~ ., data = train_data,  weights=fcbps$weights)
Yevalhat<- predict(fcbps_lm, as.data.frame(FPCAscores_eval))
fcbps_mse_lm_eval <- mean((Yeval - Yevalhat)^2)
print(fcbps_mse_lm_eval)



### npfcbps estimator
npfcbps_lm <- lm(Y ~ ., data = train_data, weights=npfcbps$w)
Yevalhat<- predict(npfcbps_lm, as.data.frame(FPCAscores_eval))
npfcbps_mse_lm_eval <- mean((Yeval - Yevalhat)^2)
print(npfcbps_mse_lm_eval)


### reg estimator 
lam_reg1 <- KRR_tune(KA * KX, Y, lam_seq_reg, criterion = "loo", alpha = 1) 
KXbar <- colMeans(KX)
Yfit <- KRR(KA * KX, gA = KA * KX, Y, lam = lam_reg1$select_lam)
Yfitadj <- KRR(KA * KX, gA = KA * matrix(KXbar, nrow = dim(KA)[2], ncol = ncol(KX), byrow = TRUE), Y, lam = lam_reg1$select_lam)
Yevalreg <- KRR(KA * KX, gA = t(KAa) * matrix(KXbar, nrow = dim(KAa)[2], ncol = ncol(KX), byrow = TRUE), Y, lam = lam_reg1$select_lam)
reg_mse <- mean((Yeval - Yevalreg)^2)
print(reg_mse)


####### FLM estimator #########
FPCscores <- as.data.frame(FPCscores)
colnames(X) <- c("X1", "X2", "X3", "X4")
train_data <- as.data.frame(cbind(Y, FPCscores, X))
summary(train_data)
lm_fit <- lm(Y ~ ., data = train_data)
TX <- mean(X %*% lm_fit$coefficients[c("X1", "X2", "X3", "X4")] )
coef <- as.numeric(lm_fit$coefficients[!names(lm_fit$coefficients) %in% c("(Intercept)","X1", "X2", "X3", "X4")])
Yevalhat <- lm_fit$coefficients[1] + as.matrix(FPCAscores_eval) %*% coef + TX
lm_eval_mse <- mean((Yeval - Yevalhat)^2)
print(lm_eval_mse)



####### FGAM estimator 
gamfit <- pfr(Y ~ af(A, k=c(6,6), m =list(c(2,3), c(2,3)), bs="ps",
presmooth="fpca.face") + s(X1, bs = "cr") + s(X2, bs = "cr") + s(X3, bs = "cr") + s(X4, bs = "cr"))
get_effect <- function(a){
Atemp <- matrix(a, nrow = nrow(X) , ncol = ncol(A), byrow = TRUE)
newdata <- list(A = Atemp, X1 = X1, X2 = X2, X3 = X3, X4 = X4)
pred <- predict(gamfit, newdata=newdata)
return(mean(pred))
}
get_mult_effect <- function(Aeval){
    Ytemp <- rep(0, nrow(Aeval))
    for (i in 1:nrow(Aeval)){
        a = Aeval[i,]
        Ytemp[i] = get_effect(a)
    }
    return(Ytemp)
}
Yevalhat <- get_mult_effect(Aeval)
gam_eval_mse <- mean((Yeval - Yevalhat)^2)
print(gam_eval_mse)


### cfb estimator 
tau_seq <- c(0, (exp(seq(log(1e-02), log(10), length.out = 10))))
lam_fcbps <- KRR_tune(KA, Y * fcbps$weights * n, lam_seq, criterion = "loo", alpha = 1)
out <- cfb_KRR(KA, KA, KA, KX, lam = lam_fcbps$select_lam, tau_seq = tau_seq, maxit = 4000, xtol_rel = 1e-8, lower =NULL)
comp_errors <- function(Y, Yfit, Yfitadj, lam, tau_idx){
        comp1 =   KRR(KA, gA = KA, Yfit * out$wlist[tau_idx,], lam = lam) - Yfitadj
        comp2 = KRR(KA, gA = KA, (Y-Yfit) * out$wlist[tau_idx,], lam = lam) 
        return(mean(comp1^2) + mean(comp2^2))
    }
test_errors <- rep(0, length(tau_seq))
for(i in 1:length(tau_seq)){
    test_errors[i] = comp_errors(Y, Yfit, Yfitadj, lam_fcbps$select_lam, tau_idx = i)
}
select_tau_idx <- which.min(test_errors)
Yevalhat <- KRR(KA, gA = t(KAa), Y * out$wlist[select_tau_idx,],lam = lam_fcbps$select_lam)
cfb_mse <- mean((Yeval - Yevalhat)^2)
print(cfb_mse)


