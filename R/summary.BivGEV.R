
################################################################################################
#################################   summary.BivGEV     #########################################
################################################################################################

summary.BivGEV <- function (object, n.sim = 1000, prob.lev = 0.05)
{
  
  testStat <- getFromNamespace("testStat", "mgcv")
  liu2 <- getFromNamespace("liu2", "mgcv")
  bs <- SE <- Vb <- est.RHOb <- V <- 1
  n <- object$n
  n.sel <- object$n.sel
  VC <- object$VC
  tau.eq1 <- object$tau.eq1
  tau.eq2 <- object$tau.eq2
  formula <- object$formula
  
  tableN <- table <- list(NULL, NULL)
  CIkt <- CInu <- NULL
  epsilon <- 1e-07
  max.p <- 0.9999999
  
  est.RHOb <- rep(NA, n.sim)
  lf <- length(object$fit$argument)
  coefficients <- object$fit$argument
  Vb <- object$Vb
  SE <- sqrt(diag(Vb))
  
  bs <- rMVN(n.sim, mean = coefficients, sigma = Vb)
  epds <- bs[, lf]
  
  ############ teta #############################
  
  est.RHOb <- teta.tr(object$VC, epds)$teta
  est.RHOb <- t(as.matrix(est.RHOb))
  CIrs <- rowQuantiles(est.RHOb, probs = c(prob.lev/2, 1 - prob.lev/2), na.rm = TRUE)
  CIrs <- t(CIrs)
  
  ############ tau #############################
  
  if (object$VC$BivD %in% c("J0", "J180"))
    est.RHOb <- ifelse(est.RHOb > 50, 50, est.RHOb)
  if (object$VC$BivD %in% c("J90", "J270"))
    est.RHOb <- -ifelse(abs(est.RHOb) > 50, 50, abs(est.RHOb))
  if (object$VC$BivD %in% c("C0", "C180", "G0", "G180"))
    est.RHOb <- ifelse(est.RHOb > 100, 100, est.RHOb)
  if (object$VC$BivD %in% c("C90", "C270", "G90", "G270"))
    est.RHOb <- -ifelse(abs(est.RHOb) > 100, 100, abs(est.RHOb))
  if (!(object$VC$BivD %in% c("AMH", "FGM")))
    tau <- BiCopPar2Tau(family = object$VC$nCa, par = est.RHOb)
  if (object$VC$BivD == "AMH")
    tau <- 1 - (2/3)/est.RHOb^2 * (est.RHOb + (1 - est.RHOb)^2 * log(1 - est.RHOb))
  if (object$VC$BivD == "FGM")
    tau <- 2/9 * est.RHOb
  CIkt <- rowQuantiles(tau, probs = c(prob.lev/2, 1 - prob.lev/2), na.rm = TRUE)
  CIkt <- t(CIkt)
  
  tau.a = object$tau.a
  
  ############ estimate ################################
  
  index <- 1:2
  ind1 <- 1:VC$X1.d2
  ind2 <- (VC$X1.d2 + 1):(VC$X1.d2 + VC$X2.d2)
  
  ind <- list(ind1 = ind1, ind2 = ind2)
  for (i in index) {
    estimate <- round(object$fit$argument[ind[[i]]],5)
    se <- SE[ind[[i]]]
    ratio <- round((estimate/se),5)
    pv <- round(2 * pnorm(abs(ratio), lower.tail = FALSE),5)
    table[[i]] <- cbind(estimate, se, ratio, pv)
    dimnames(table[[i]])[[2]] <- c("Estimate", "Std. Error",
                                   "z value", "Pr(>|z|)")
  }
  
  res <- list(tableP1 = table[[1]], tableP2 = table[[2]],
              formula = formula, n = n, n.sel=n.sel, CIrs = CIrs, CIkt = CIkt, tau = object$tau, tau.a = object$tau.a,
              theta.a = object$theta.a, theta = object$theta, t.edf = object$t.edf, BivD = object$BivD, Model = object$Model,
              tau.eq1 = tau.eq1, tau.eq2 = tau.eq2)
  
  class(res) <- "summary.BivGEV"
  res
}

