

#######################################################################################################
######################################     BivGEV       ###############################################
#######################################################################################################

BivGEV <- function(formula, data = list(), 
    Model = "BivGEV", BivD = "N", rinit = 1, rmax = 100, iterlim = 100, tolsp = 1e-07,
    parscale, tau.eq1 = -0.25, tau.eq2= -0.25, gc.m=FALSE, max.pr = 0.999999)
{

i.rho <- y1.y2 <- y1.cy2 <- cy1.y2 <- cy1.cy2 <- cy <- cy1 <- inde <- n.sel <- NULL  #q.mag <- sp
end <- 0  #l.sp1 <- l.sp2 <- l.sp3
X2s <- NULL
i.rho <- 0
mod1 <- mod2 <- NULL
weights <- NULL
#ngc <- 2
opc <- c("N", "C0", "C90", "C180", "C270", "J0", "J90", "J180",
        "J270", "G0", "G90", "G180", "G270", "F", "AMH", "FGM")
scc <- c("C0", "C180", "J0", "J180", "G0", "G180")
sccn <- c("C90", "C270", "J90", "J270", "G90", "G270")

 l.flist <- length(formula)
if (!is.list(formula))
    stop("You must specify a list of equations.")

cl <- match.call()
mf <- match.call(expand.dots = FALSE)
cha2 <- as.character(formula[[2]][2])
ig <- interpret.gam(formula)
if (l.flist == 2) {
        e1 <- all.vars(as.formula(formula[[1]]))[1]
        e1 <- c(e1, ig[[1]]$pred.names)
        e2 <- all.vars(as.formula(formula[[2]]))[1]
        e2 <- c(e2, ig[[2]]$pred.names)
    pred.n <- union(e1, c(e2, cha2))
  }

fake.formula <- paste(e1[1], "~", paste(pred.n, collapse = " + "))
environment(fake.formula) <- environment(formula[[1]])
mf$formula <- fake.formula
mf$Model <- mf$BivD <-  mf$rinit <- mf$rmax <- mf$iterlim <- mf$tolsp <- mf$gc.m <- mf$tau.eq1 <- mf$tau.eq2 <- mf$parscale <- NULL  #<- mf$dof  
mf$drop.unused.levels <- TRUE


if(Model=="SampleSelGEV") mf$na.action <- na.pass   ####### SAMPLE SELECTION

mf[[1]] <- as.name("model.frame")
data    <- eval(mf, parent.frame())

if(gc.m == TRUE) gc()

if(Model=="SampleSelGEV"){                        ####### SAMPLE SELECTION
  data[is.na(data[, e1[1]]), e1[1]] <- 0
  indS <- data[, e1[1]]
  indS[is.na(indS)] <- 0
  indS <- as.logical(indS)
  data[indS == FALSE, e2[1]] <- 0
  data <- na.omit(data)
}


if(!("(weights)" %in% names(data))) {weights <- rep(1,dim(data)[1]) 
data$weights <- weights
names(data)[length(names(data))] <- "(weights)"} else weights <- data[,"(weights)"] 


formula.eq1 <- formula[[1]]
formula.eq2 <- formula[[2]]

if (Model == "BivGEV") {
    if (e1[1] %in% e2[-1])
    end <- 1
    if (e2[1] %in% e1[-1])
      end <- 2
  }

ct  <- data.frame(c(opc), c(1:14, 55, 56))
cta <- data.frame(c(opc), c(1, 3, 23, 13, 33, 6, 26, 16, 36, 4, 24, 14, 34, 5, 55, 56))
nC  <- ct[which(ct[, 1] == BivD), 2]
nCa <- cta[which(cta[, 1] == BivD), 2]




##############################################################  
# Equation 1                                                 #
############################################################## 


mod1 <- eval(substitute(bgeva(formula.eq1, tau = tau.eq1, data = data, Hes=T)))

X1    <- mod1$X
X1.d2 <- dim(X1)[2]
y1    <- mod1$gam.fit$y
n     <- length(y1)
inde  <- rep(TRUE, n)



##############################################################  
# Equation 2                                                 #
############################################################## 

if (Model == "BivGEV") {
mod2 <- eval(substitute(bgeva(formula.eq2, tau = tau.eq2, data = data,  Hes=T)))
X2    <- mod2$X
X2.d2 <- dim(X2)[2]
y2 <- mod2$gam.fit$y
#        if (Model == "BivGEV") {
            y1.y2 <- y1 * y2
            y1.cy2 <- y1 * (1 - y2)
            cy1.y2 <- (1 - y1) * y2
            cy1.cy2 <- (1 - y1) * (1 - y2)
#        }
}

if (Model == "BivGEV") {
beta1 <- c(coef(mod1))
beta2 <- c(coef(mod2))

eta1 <- X1%*%beta1
eta2 <- X2%*%beta2

p1<- exp(-( 1 + tau.eq1*(X1%*%beta1))^(-1/tau.eq1))
p2<- exp(-( 1 + tau.eq2*(X2%*%beta2))^(-1/tau.eq2))

y1m    <- mod1$gam.fit$y
y1star <- ifelse(y1m==1,1,-1)
res.dev1 <-  y1star * (sqrt(-2*(y1m*(log(p1)) + (1-y1m)*(log(1-p1)))))

y2m    <- mod2$gam.fit$y
y2star <- ifelse(y1m==1,1,-1)
res.dev2 <-  y2star * (sqrt(-2*(y2m*(log(p2)) + (1-y2m)*(log(1-p2)))))


        res1 <- res.dev1
        res2 <- res.dev2
        ass.s <- cor(res1, res2, method = "kendall")
        ss <- sign(ass.s)
        ass.s <- ss*ifelse(abs(ass.s) > 0.9, 0.9, abs(ass.s))      
}



##############################################################  
# Equation 2 Sample Selection                                #
############################################################## 

if (Model == "SampleSelGEV") {
    inde <- as.logical(y1)
    mod2 <- eval(substitute(bgeva(formula.eq2, tau = tau.eq2, data = data[inde,],  Hes=T)))
                 

    ##################################
    ######### definisco X2s ##########
    ##################################
    
    mod2.s <- eval(substitute(bgeva(eq2, tau = tau.eq2, data = data,  Hes=T)))
    X2s <- mod2.s$X    
    
#####################
X2     <- mod2$X
X2.d2  <- dim(X2)[2]
y2     <- mod2$gam.fit$y
n.sel  <- sum(as.numeric(inde))
cy1    <- (1-y1)
y1.y2  <- y1[inde]*y2
y1.cy2 <- y1[inde]*(1-y2)
######################
epsilon <- 1e-07
max.p <- 0.9999999

formula.eq2imr <- update.formula(formula.eq2, ~. + imr)
beta1 <- c(coef(mod1))
eta1 <- X1%*%beta1
p.s1 <- exp(-( 1 + tau.eq1*(X1%*%beta1))^(-1/tau.eq1))
p.s1 <- ifelse(p.s1 < epsilon, epsilon, p.s1)
p.s1 <- ifelse(p.s1 > max.p, max.p, p.s1)
#p.s1 <- ifelse(is.na(p.s1), quantile(p.s1, c(0.5), na.rm=T), p.s1)

data$imr <- dnorm(p.s1)/pnorm(p.s1)
data$imr <- ifelse(is.na(data$imr), quantile(data$imr, c(0.5), na.rm=T), data$imr)
imr <- data$imr
      
mod2.1 <- eval(substitute(bgeva(formula.eq2imr, tau = tau.eq2, data = data[inde,],Hes=T)))
      
pimr <- which(names(coef(mod2.1))=="imr")
c.mod2 <- coef(mod2.1)[-pimr]

beta2s <- c(coef(mod2.1))
eta2.s <- mod2.1$X%*%beta2s
p.s2 <- exp(-( 1 + tau.eq2*(eta2.s))^(-1/tau.eq2))
p.s2 <- ifelse(p.s2 < epsilon, epsilon, p.s2)
p.s2 <- ifelse(p.s2 > max.p, max.p, p.s2)
#p.s2 <- ifelse(is.na(p.s2), quantile(p.s2, c(0.5), na.rm=T), p.s2)

y2.1    <- mod2.1$gam.fit$y
ystar <- ifelse(y2.1==1,1,-1)
rr2 <-  ystar * (sqrt(-2*(y2.1*(log(p.s2)) + (1-y2.1)*(log(1-p.s2)))))
#rr2 <- ifelse(is.na(rr2) | rr2==-Inf | rr2==Inf, quantile(rr2, c(0.5), na.rm=T), rr2)

sia <- sqrt( (mean(rr2)^2) + mean(imr[inde]*(imr[inde]+p.s1[inde]))*mod2.1$coef["imr"]^2)[[1]]
co <- (mod2.1$coef["imr"]/sia)[[1]]
ass.s <- co
ss <- sign(ass.s)
ass.s <- ss*ifelse(abs(ass.s) > 0.2, 0.2, abs(ass.s))
      
 }

        if (BivD %in% scc)
        ass.s <- abs(ass.s)
        if (BivD %in% sccn)
        ass.s <- -abs(ass.s)
        if (!(BivD %in% c("AMH", "FGM")))
        i.rho <- BiCopTau2Par(family = nCa, tau = ass.s, check.taus=F)
        if (BivD %in% c("AMH", "FGM"))
        i.rho <- BiCopTau2Par(family = 1, tau = ass.s, check.taus=F)
        if (BivD %in% c("N", "AMH", "FGM"))
        i.rho <- atanh(i.rho)
        if (BivD == "F")
        i.rho <- ifelse(abs(i.rho) < 1e-07, 1e-07, i.rho)
        if (!(BivD %in% c("N", "AMH", "FGM", "F")))
        i.rho <- abs(i.rho)
        if (BivD %in% c("C0", "C180", "C90", "C270"))
        i.rho <- log(i.rho)
        if (BivD %in% c("J0", "J180", "G0", "G180", "J90", "J270", "G90", "G270"))
        i.rho <- log(i.rho - 1)
      

names(i.rho) <- "theta.star"
if(Model == "BivGEV") start.v <- c(coef(mod1), coef(mod2), i.rho)
if(Model == "SampleSelGEV") start.v <- c(coef(mod1), c.mod2, i.rho)


if (missing(parscale))
parscale <- 1

respvec <- list(y1 = y1, y2 = y2, y1.y2 = y1.y2, y1.cy2 = y1.cy2,
                  cy1.y2 = cy1.y2, cy1.cy2 = cy1.cy2, cy1 = cy1, cy = cy,
                  univ = 0)
respvec1 <- respvec
respvec1$univ <- 1

VC <- list(X1 = X1, X2 = X2,
           X1.d2 = X1.d2,
           X2.d2 = X2.d2,
           X2s = X2s,
           weights = weights,
           #Hes = Hes,
           nCa = nCa,
           end = end,
           BivD = BivD,
           nC = nC,
           n = n,
           inde = inde,
           parscale = parscale,
           tau.eq1=tau.eq1,
           tau.eq2=tau.eq2,
           nCa = nCa,
           gc.m = gc.m,
           Model = Model)

if (Model == "BivGEV")   func.opt <- BivGEVOptim
if (Model == "SampleSelGEV") func.opt <- BivGEVOptimSS



BivGEVFit   <- BivGEV.fit(func.opt = func.opt, start.v = start.v, rinit = rinit, rmax = rmax, 
                          iterlim = 100, respvec = respvec, VC = VC)

BivGEVFit.p <- BivGEVFit.post(BivGEVFit = BivGEVFit, Model = Model, VC = VC, mod1 = mod1, mod2 = mod2)
BivGEVFit   <- BivGEVFit.p$BivGEVFit

if(gc.m == TRUE) 
  gc()


rm(respvec1)

##### check da qui in poi
e.v <- min(eigen(BivGEVFit$fit$hessian, symmetric = TRUE, only.values = TRUE)$values)
e.v
gradi <- round(max(abs(BivGEVFit$fit$gradient)), 2)
gradi
me1 <- "Largest absolute gradient value is not close to 0."
me2 <- "Information matrix is not positive definite."
me3 <- "Read the WARNINGS section in ?BivGev."
a <- if (gradi > 0.1 && e.v <= 0) {
    warning(me1, call. = FALSE)
    warning(paste(me2, "\n", me3), call. = FALSE)
    }
b <- if (gradi > 0.1 && e.v > 0)
    warning(paste(me1, "\n", me3), call. = FALSE)
c <- if (gradi < 0.1 && e.v <= 0)
    warning(paste(me2, "\n", me3), call. = FALSE)

dataset <- data
rm(data)

L <- list(fit = BivGEVFit$fit, 
          coefficients = BivGEVFit$fit$argument, 
          p00 = BivGEVFit$fit$p00, 
          p01 = BivGEVFit$fit$p01, 
          p11 = BivGEVFit$fit$p11, 
          p10 = BivGEVFit$fit$p10, 
          p1 = BivGEVFit$fit$p1, 
          p2 = BivGEVFit$fit$p2, 
          tau.eq1 = tau.eq1, 
          tau.eq2 = tau.eq2, 
          t.edf = BivGEVFit.p$t.edf, 
          bgev1 = mod1, 
          bgev2 = mod2, 
          iter.if = BivGEVFit$iter.if,
          bs.mgfit = BivGEVFit$bs.mgfit, 
          post.f = BivGEVFit$post.f, 
          magpp = BivGEVFit$magpp,
          He = BivGEVFit.p$He, 
          logLik = BivGEVFit.p$logLik, 
          Vb = BivGEVFit.p$Vb,
          HeSh = BivGEVFit.p$HeSh, 
          He = BivGEVFit.p$He, 
          F = BivGEVFit.p$F, 
          theta = BivGEVFit.p$theta, 
          theta.a = BivGEVFit.p$theta.a, 
          #dof = BivGEVFit.p$dof,
          tau = BivGEVFit.p$tau, 
          tau.a = BivGEVFit.p$tau.a, 
          R = BivGEVFit.p$R, 
          Ve = BivGEVFit.p$Ve,
          formula = formula, 
          n = n, 
          n.sel = n.sel, 
          X1 = X1, 
          X2 = X2, 
          X1.d2 = X1.d2, 
          X2.d2 = X2.d2, 
          X2s = X2s,
          nC = nC, 
          nCa = nCa, 
          BivD = BivD, 
          VC = VC, 
          respvec = respvec, 
          l.flist = l.flist, 
          dataset = dataset,
          Model=Model,
          y1 = y1,
          y2 = y2)
class(L) <- "BivGEV"
L
}


