\name{BivGEV}
\alias{BivGEV}

\title{
Bivariate Generalised Extreme Value Model
}
\description{
BivGEV function is a wrapper function from \code{GJRM} R package. It fits joint models with imbalanced binary margins.

}
\usage{
BivGEV(formula, data = list(), Model = "BivGEV", 
BivD = "N", rinit = 1, rmax = 100, iterlim = 100, tolsp = 1e-07, 
parscale, tau.eq1 = -0.25, tau.eq2 = -0.25, gc.m=FALSE, 
min.dn = 1e-40, min.pr = 1e-16, max.pr = 0.999999)
}

\arguments{
  \item{formula}{
The formula will be a list of two equations.
}
  \item{data}{
The data frame includes the variables of the model. 
}
  \item{Model}{
Possible choice are "BivGEV" (Bivariate Generalised Extreme Value) or "SampleSelGEV" (Sample Selection Generalised Extreme Value).
}
  \item{BivD}{
The possible typologies of bivariate error distribution employed are "N" (Gaussian), "C0" (Clayton), "J0" (Joe), "G0" (Gumbel), or "F" (Frank). 
}

  \item{rinit}{
Starting the trust-region radius.
}
  \item{rmax}{
Maximum allowed trust region radius.
}
  \item{iterlim}{
A positive integer specifying the maximum number of loops.
}
  \item{tolsp}{
Tolerance for convergence of the algorithm.
}

  \item{parscale}{
Parscale is a scalar parameter. If the option is missing then any re-scaling is done during the model fitting process.
}

  \item{tau.eq1}{
The shape parameter of the GEV distribution for the first equation. 
}
  \item{tau.eq2}{
Shape parameter of the GEV distribution for the second equation. 
}

 \item{gc.m}{
This option is useful when considering large data samples. If TRUE then the garbage collector is taken into account. 
}

 \item{min.dn}{The minimum threshold of densities during the model estimation process.
}

 \item{min.pr}{The minimum threshold of probabilities during the model estimation process.
}

 \item{max.pr}{The maximum threshold of probabilities during the model estimation process.
}

}
\details{
The Bivariate Generalized Extreme Value model incorporates the quantile function of the Generalized Extreme Value (GEV) random variable as link function for each equation of the model. It considers the dependence between the two equations throught the use of a Copula function.

The shape parameter \eqn{\tau} of the GEV distribution controls the tail behaviour. For different values of the parameter \eqn{\tau}, three families of distributions are defined:
\itemize{ 
\item for \eqn{\tau} --> 0, the GEV distribution is known as the Gumbel class; 
\item for \eqn{\tau} > 0, the GEV distribution is known as Frechet;
\item for \eqn{\tau} < 0, the GEV distribution is known as Weibull.
}
For example, if the sample is imbalanced with a low percentage of Y = 1, as in credit scoring models, the best link function is obtained for \eqn{\tau} < 0.

More details can be found in Calabrese, Osmetti and Zanin (2019).
}
\value{
The function returns an object of class BivGEV.
}
\references{
Calabrese R., Osmetti S.A., Zanin L. (2019). A joint scoring model for peer-to-peer and traditional lending: a bivariate model with copula dependence. \emph{Journal of the Royal Statistical Society Serie A}, 182, 1163-1188 

Calabrese R., Osmetti S.A., Zanin L. (2020). Sample selection bias in peer-to-peer lending market. 

Marra G. and Radice R. (2019). \emph{GJRM: Generalised Joint Regression Modelling}. R package Version: 0.2
}

\author{
Raffaella Calabrese - \email{raffaella.calabrese@ed.ac.uk};
Silvia Angela Osmetti - \email{silvia.osmetti@unicatt.it}; 
Maintainer: Luca Zanin - \email{luca.zanin@studio.unibo.it};
}

\note{
Some issues of convergence may occur because some combinations of \emph{tau.eq1} and \emph{tau.eq2} might not be supported by data. Convergence can be controlled using \code{\link{diagnostics}} function.
}

\seealso{
\code{\link{diagnostics}}, \code{\link{predict.BivGEV}}, \code{\link{summary.BivGEV}}, \code{\link{train.test}}
}

\examples{

### library(devtools)
### Install the BivGEV package from the Github repository ###
### install_github("BivGEV/BivGEV")
 
library(BivGEV)
library(evd)
library(VineCopula)
library(mgcv)

#################################################################
########### Bivariate Generalised Extreme Value Model ########### 
#################################################################

#################################
#### Data Generating Process ####
#################################

set.seed(0)

set.theta <- BiCopTau2Par(14, 0.5, check.taus = TRUE)

n      <- 10000      # Sample size
family <- 14        # Survival Gumbel copula
theta  <- set.theta # Copula parameter

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

##################################################################
#### Example 1:                                               ####
#### Estimation of a BivGEV with Gumbel Copula                ####
##################################################################

#### Fixed tau.eq1 and tau.eq2
tau1 <- -0.30
tau2 <- -0.30

eq1 <- y1 ~ x 
eq2 <- y2 ~ x

out1 <- BivGEV(list(eq1,eq2), data = dataSim, BivD = "G0", 
       Model = "BivGEV", tau.eq1 = tau1, tau.eq2 = tau2)
summary(out1)
diagnostics(out1)






########################################################################
########### Sample selection Generalised Extreme Value Model ########### 
########################################################################

#################################
#### Data Generating Process ####
#################################

set.seed(10)

set.theta <- BiCopTau2Par(14, 0.5, check.taus = TRUE)

n      <- 50000      # Sample size
family <- 14         # Survival Gumbel copula
theta  <- set.theta  # Copula parameter

U  <- BiCopSim(n, family, theta)
ER <- cbind(qgev(1-U[,1], shape=-0.30),qgev(1-U[,2],shape=-0.35))

sigma       <- matrix(0.5,3, 3) 
diag(sigma) <- 1
covariance  <- rmvn(n, rep(0,3), sigma)
covariance  <- pnorm(covariance)
x1 <- covariance[,1]
x2 <- covariance[,2]
x3 <- covariance[,3]

eta1 <-  -1.2 + 1.7*x1 - 0.9*x2 + 1.3*x3 
eta2 <-  -1.7 + 1.3*x1 - 0.1*x2 

ysel <-  ifelse(eta1 - ER[, 1] > 0,1,0)
y    <-  ifelse(eta2 - ER[, 2] > 0,1,0)
yobs <-  y*(ysel > 0)

dataSim <- data.frame(ysel,yobs,x1,x2,x3)


############################################################################
#### Example 2:                                                         ####
#### Estimation of a BivGEV Sample selection model with Gumbel Copula   ####
############################################################################

#### Fixed tau.eq1 and tau.eq2
tau1 <- -0.30
tau2 <- -0.35

eq1 <- ysel ~ x1 + x2 + x3 
eq2 <- yobs ~ x1 + x2

out2 <- BivGEV(list(eq1,eq2), data = dataSim, BivD = "G0", 
       Model = "SampleSelGEV", tau.eq1 = tau1, tau.eq2 = tau2)
summary(out2)

  }


\keyword{Bivariate Generalised Extreme Value Model}
