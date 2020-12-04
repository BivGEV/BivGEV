

##########################################################################################
##########################           predict.BivGEV         ##############################
##########################################################################################

predict.BivGEV <- function (object, newdata, type = "joint", ...)
{
  
  if(type == "naive"){
    ss.pred1 <- names(object$bgev1$coefficients)
    tau.eq1  <- object$tau.eq1
    formula1 <- as.formula(object$formula[[1]])
    X1       <- model.matrix(formula1, data=newdata)
    eta1     <- X1 %*% object$bgev1$coefficients
    p1.naive <- as.numeric(exp(-(1 + tau.eq1 * (eta1))^(-1/tau.eq1)))
    
    ss.pred2 <- names(object$bgev2$coefficients)
    tau.eq2  <- object$tau.eq2
    formula2 <- as.formula(object$formula[[2]])
    X2       <- model.matrix(formula2, data=newdata)
    eta2     <- X2 %*% object$bgev2$coefficients
    p2.naive <- as.numeric(exp(-(1 + tau.eq2 * (eta2))^(-1/tau.eq2)))
  }
  
  
  if(type == "joint"){
    
    ss.pred1 <- names(object$bgev1$coefficients)
    ind1     <- 1:object$X1.d2
    tau.eq1  <- object$tau.eq1
    formula1 <- as.formula(object$formula[[1]])
    X1       <- model.matrix(formula1, data=newdata)
    eta1     <- X1 %*% object$coefficients[ind1]
    p1       <- as.numeric(exp(-(1 + tau.eq1 * (eta1))^(-1/tau.eq1)))
    
    ss.pred2 <- names(object$bgev2$coefficients)
    ind2     <- (object$X1.d2 + 1):(object$X1.d2 + object$X2.d2)
    tau.eq2  <- object$tau.eq2
    formula2 <- as.formula(object$formula[[2]])
    X2       <- model.matrix(formula2, data=newdata)
    eta2     <- X2 %*% object$coefficients[ind2]
    p2       <- as.numeric(exp(-(1 + tau.eq2 * (eta2))^(-1/tau.eq2)))
    theta    <- object$theta
    epsilon  <- 1e-07
    
    p11 <- BiCDF(p1, p2, object$nC, theta)
    p10 <- pmax(p1 - p11, epsilon)
    p01 <- pmax(p2 - p11, epsilon)
    p00 <- pmax(1 - p11 - p10 - p01, epsilon)
    
  }
  
  if(type == "naive"){
    result <- data.frame(p1.naive, p2.naive)
    return(result)
  }
  
  if(type == "joint"){
    result <- data.frame(p1, p2, p11, p10, p01, p00)
    return(result)
  }
}


