


################################################################################################
##############################     BivGEV.fit.post      ########################################
################################################################################################

BivGEVFit.post <- function (BivGEVFit, Model, VC, mod1, mod2) 
{

  Ve <- X2s <- theta <- theta.a <- p1n <- p2n <- R <- NULL
  He <- BivGEVFit$fit$hessian
  logLik <- -BivGEVFit$fit$l
  epsilon <- 1e-07
  max.p <- 0.9999999
  He.eig <- eigen(He, symmetric = TRUE)
  
  if (min(He.eig$values) < sqrt(.Machine$double.eps) && sign(min(sign(He.eig$values))) == -1)
    He.eig$values <- abs(He.eig$values)
  
  if (min(He.eig$values) < sqrt(.Machine$double.eps)) {
    pep <- which(He.eig$values < sqrt(.Machine$double.eps))
    He.eig$values[pep] <- epsilon
  }
  
  Vb <- He.eig$vectors %*% tcrossprod(diag(1/He.eig$values), He.eig$vectors)
  Vb <- (Vb + t(Vb))/2
  
  HeSh <- He
  Ve <- Vb
  F <- diag(rep(1, dim(Vb)[1]))
  R <- BivGEVFit$bs.mgfit$R
  
  t.edf <- sum(diag(F))
  
  dimnames(BivGEVFit$fit$hessian)[[1]] <- dimnames(BivGEVFit$fit$hessian)[[2]] <- dimnames(Vb)[[1]] <- dimnames(Vb)[[2]] <- dimnames(HeSh)[[1]] <- dimnames(HeSh)[[2]] <- dimnames(F)[[1]] <- dimnames(F)[[2]] <- dimnames(He)[[1]] <- dimnames(He)[[2]] <- names(BivGEVFit$fit$argument)
  

  if (Model == "BivGEV") {
  dep <- BivGEVFit$fit$etad
  names(dep) <- "theta"
  theta <- teta.tr(VC, dep)$teta
  theta.a <- mean(theta)
  }

  
  if (Model == "BivGEVss") {
    BivGEVFit$fit$eta2 <- VC$X2s %*% BivGEVFit$fit$argument[(VC$X1.d2 + 1):(VC$X1.d2 + VC$X2.d2)]
    
    tau.eq1 <- VC$tau.eq1
    beta1.1 <- c(coef(mod1))
    eta1.0 <- VC$X1 %*%beta1.1 
    p1n <- exp(-(1 + tau.eq1 * (eta1.0))^(-1/tau.eq1)) ## naive
    p1n <- ifelse(p1n < epsilon, epsilon, p1n)
    p1n <- ifelse(p1n > max.p, max.p, p1n)
  
    tau.eq2 <- VC$tau.eq2
    beta2.2 <- c(coef(mod2))
    eta2.0 <- VC$X2s %*%beta2.2
    p2n <- exp(-(1 + tau.eq2 * (eta2.0))^(-1/tau.eq2)) ## naive
    p2n <- ifelse(p2n < epsilon, epsilon, p2n)
    p2n <- ifelse(p2n > max.p, max.p, p2n)  
  
      dep <- BivGEVFit$fit$etad
      names(dep) <- "theta"
      theta <- teta.tr(VC, dep)$teta
      theta.a <- mean(theta)
  }
  
  
  
  if (Model == "BivGEVss") {
  
    p1 <- exp(-(1 + tau.eq1 * (BivGEVFit$fit$eta1))^(-1/tau.eq1))
    p1 <- ifelse(p1 < epsilon, epsilon, p1)
    p1 <- ifelse(p1 > max.p, max.p, p1)
    
    p2 <- exp(-(1 + tau.eq2 * (BivGEVFit$fit$eta2))^(-1/tau.eq2)) # eta2 modified above
    p2 <- ifelse(p2 < epsilon, epsilon, p2)
    p2 <- ifelse(p2 > max.p, max.p, p2)  

    p11 <- BiCDF(p1, p2, VC$nC, theta) #, VC$dof)
    BivGEVFit$fit$p10 <- p1 - p11
    BivGEVFit$fit$p11 <- p11
    BivGEVFit$fit$p00 <- (1 - p2) - (p1 - p11)
    BivGEVFit$fit$p01 <- p2 - p11
    BivGEVFit$fit$p1 <- p1
    BivGEVFit$fit$p2 <- p2
    
  }
  
 
  if (VC$BivD %in% c("J0", "J180", "J90", "J270"))
    theta <- ifelse(abs(theta) > 50, 50, abs(theta))
  if (VC$BivD %in% c("J90", "J270"))
    theta <- -theta
  if (VC$BivD %in% c("C0", "C180", "G0", "G180", "C90", "C270",
                     "G90", "G270"))
    theta <- ifelse(abs(theta) > 100, 100, abs(theta))
  if (VC$BivD %in% c("C90", "C270", "G90", "G270"))
    theta <- -theta
  if (!(VC$BivD %in% c("AMH", "FGM")))
    tau <- BiCopPar2Tau(family = VC$nCa, par = theta)
  if (VC$BivD == "AMH")
    tau <- 1 - (2/3)/theta^2 * (theta + (1 - theta)^2 * log(1 - theta))
  if (VC$BivD == "FGM")
    tau <- 2/9 * theta
  tau.a <- mean(tau)
  if (length(tau) == 1)
    tau.a <- tau
  
  if(Model == "BivGEV"){
  p00 <- BivGEVFit$fit$p00
  p01 <- BivGEVFit$fit$p01
  p11 <- BivGEVFit$fit$p11
  p10 <- BivGEVFit$fit$p10
  p1 <-  BivGEVFit$fit$p1
  p2 <-  BivGEVFit$fit$p2
}
  
  if (length(theta) == 1)
  theta.a <- theta  
  
  if (VC$gc.m == TRUE) 
    gc()
  
  list(BivGEVFit = BivGEVFit, He = He, logLik = logLik, Vb = Vb, 
       HeSh = HeSh, F = F, theta = theta, t.edf = t.edf, 
       theta.a = theta.a, tau = tau, tau.a = tau.a, p1n = p1n, p2n = p2n, Ve = Ve) #, dof = VC$dof, dof.a = VC$dof)

}


