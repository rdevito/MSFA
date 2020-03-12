## ----help, echo = TRUE, results = TRUE, tidy = TRUE----------------------
library(MSFA)
data(data_immune)
help(data_immune)

## ----posterior, messages = FALSE-----------------------------------------
set.seed(1971)
out10_1010 <- sp_msfa(data_immune,  k = 10,  j_s = c(10, 10), trace = FALSE)

## ----Phi-----------------------------------------------------------------
p <- ncol(data_immune[[1]])
nrun <- dim(out10_1010$Phi)[3]
SigmaPhi <-  SigmaLambda1 <- SigmaLambda2 <- array(0, dim=c(p, p, nrun))
for(j in 1:nrun)
{
  SigmaPhi[,,j] <- tcrossprod(out10_1010$Phi[,,j]) 
  SigmaLambda1[,,j] <- tcrossprod(out10_1010$Lambda[[1]][,,j]) 
  SigmaLambda2[,,j] <- tcrossprod(out10_1010$Lambda[[2]][,,j]) 
}
SigmaPhi <- apply(SigmaPhi, c(1, 2), median)
SigmaLambda1 <- apply(SigmaLambda1, c(1, 2), median)
SigmaLambda2 <- apply(SigmaLambda2, c(1, 2), median)
Phi <- apply(out10_1010$Phi, c(1, 2), median)

## ----SigmaPhi, fig.width=4.5, fig.height=4.5-----------------------------
plot(sp_eigen(SigmaPhi), pch = 16)
abline(h = 0.05, col = 2)

## ----SigmaLambda, fig.width=6.5, fig.height=4.5--------------------------
par(mfrow=c(1, 2))
plot(sp_eigen(SigmaLambda1), pch=16)
abline(h = 0.05, col = 2)
plot(sp_eigen(SigmaLambda2), pch=16)
abline(h = 0.05, col = 2)

## ----OP------------------------------------------------------------------
Phi_OP10 <- sp_OP(out10_1010$Phi[,1:5,], itermax = 10)

