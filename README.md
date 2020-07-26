
<img src="https://github.com/BivGEV/BivGEV/blob/master/BivGEV_logo.png" width="200">

The `BivGEV` package allows the estimation of the following models:
- Bivariate Generalised Extreme Value models
- Bivariate Generalised Extreme Value models in the presence of sample selection

The package allow the estimation of models considering differents copula, including Gaussian, Clayton, Gumbel, Joe, and Frank

## Installation

You can install the released version of BivGEV from [GitHub](https://github.com/BivGEV/BivGEV) with:

``` r
library(devtools)
install_github("BivGEV/BivGEV")  

library(BivGEV)
library(evd)
library(VineCopula)
library(mgcv)
 ```
 
##  Bivariate Generalised Extreme Value Model 

Data Generating Process
 ``` r
set.seed(0)

set.theta <- BiCopTau2Par(14, 0.5, check.taus = TRUE)

n      <- 10000      # Sample size
family <- 14         # Survival Gumbel copula
theta  <- set.theta  # Copula parameter

U  <- BiCopSim(n, family, theta)
ER <- cbind(qgev(1-U[,1], shape=-0.3),qgev(1-U[,2],shape=-0.3))

x <- rnorm(n, 1.2, 0.25)

par.eq1 <- c(-1.7, 0.4)
par.eq2 <- c(-2.2, 0.8)

eta1 <- par.eq1[1] + par.eq1[2]*x 
eta2 <- par.eq2[1] + par.eq2[2]*x

y1 <- ifelse(eta1 - ER[, 1] > 0, 1,0) 
y2 <- ifelse(eta2 - ER[, 2] > 0, 1,0)

dataSim <- data.frame(y1,y2,x)
```

Setting the tau parameter 
 ``` r
tau1 <- -0.30 # setting tau parameter for eq. 1
tau2 <- -0.30 # setting tau parameter for eq. 2
```

Model's estimation with a Survival Gumbel copula
 ``` r
eq1 <- y1 ~ x 
eq2 <- y2 ~ x
out1 <- BivGEV(list(eq1,eq2), data = dataSim, BivD = "G0", 
       Model = "BivGEV", tau.eq1 = tau1, tau.eq2 = tau2)     
summary(out1)
diagnostics(out1)
```


##  Bivariate Generalised Extreme Value model in the presence of sample selection


