\name{BivGEV}
\alias{BivGEV}

\title{
Bivariate Generalised Extreme Value Model
}
\description{
BivGEV funcion is a wrapper from \code{gjrm()} in \code{GJRM}. It fits joint models with imbalanced binary margins.

}
\usage{
BivGEV(formula, data = list(), Model = "BivGEV", 
BivD = "N", rinit = 1, rmax = 100, iterlim = 100, tolsp = 1e-07, 
parscale, tau.eq1 = -0.25, tau.eq2 = -0.25, gc.m=FALSE)
}

\arguments{
  \item{formula}{
The formula will be a list of two equations.
}
  \item{data}{
The data frame that includes the variables of the model. 
}
  \item{Model}{
Possible choice are "BivGEV" (Bivariate Generalised Extreme Value) or "BivGEVss" (Bivariate Generalised Extreme Value for Sample Selection).
}
  \item{BivD}{
The possible typologies of bivariate error distribution employed are "N" (Gaussian), "C0" (Clayton), "J0" (Joe), "G0" (Gumbel), or "F" (Frank). 
}

  \item{rinit}{
Starting trust region radius.
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
Shape parameter of the GEV distribution for the first equation. 
}
  \item{tau.eq2}{
Shape parameter of the GEV distribution for the second equation. 
}

 \item{gc.m}{
This option is useful when considering large data samples. If TRUE then the garbage collector is taken into account. 
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
Calabrese R., Osmetti S.A., Zanin L. (2019). A joint scoring model for peer-to-peer and traditional lending: a bivariate model with copula dependence. \emph{Submitted to the Journal of the Royal Statistical Society Serie A} 

%Marra G. and Radice R. (2019). \emph{GJRM: Generalised Joint Regression Modelling}. R package Version: 0.2
}

\author{
Raffaella Calabrese - \email{raffaella.calabrese@ed.ac.uk};
Silvia Angela Osmetti - \email{silvia.osmetti@unicatt.it}; 
Luca Zanin - \email{luca.zanin@studio.unibo.it};
}

\note{
Some issues of convergence may occur because some combinations of \emph{tau.eq1} and \emph{tau.eq2} might not be supported by data. Convergence can be controlled using \code{\link{diagnostics}} function.
}

\seealso{
\code{\link{diagnostics}}, \code{\link{predict.BivGEV}}, \code{\link{summary.BivGEV}}, \code{\link{train.test}}
}

\examples{

library(BivGEV)

#################################################################
########### Bivariate Generalised Extreme Value Model ########### 
#################################################################

#################################
#### Data Generating Process ####
#################################

set.seed(0)

n      <- 5000  # Sample size
family <- 1     # Gaussian copula
lambda <- 0.70  # Parameter of the copula

U  <- BiCopSim(n, family, lambda)
ER <- cbind(qnorm(U[,1]),qnorm(U[,2]))

x <- rnorm(n, 1.2, 0.25)

par.eq1 <- c(-1.8, 0.4)
par.eq2 <- c(-2.3, 0.8)

eta1 <- par.eq1[1] + par.eq1[2]*x 
eta2 <- par.eq2[1] + par.eq2[2]*x

y1 <- ifelse(eta1 + ER[, 1] > 0, 1,0) 
y2 <- ifelse(eta2 + ER[, 2] > 0, 1,0)

dataSim <- data.frame(y1,y2,x)

summary(dataSim$y1)
summary(dataSim$y2)


######################################
#### Model estimation: Example 1  ####
######################################

#### Fixed tau.eq1 and tau.eq2
tau.eq1 <- -0.325
tau.eq2 <- -0.325

eq1 <- y1 ~ x 
eq2 <- y2 ~ x

mod <- BivGEV(list(eq1,eq2), data = dataSim, BivD = "N", 
       Model = "BivGEV", tau.eq1 = tau.eq1, tau.eq2 = tau.eq2)
summary(mod)
diagnostics(mod)

#### predict marginals and Joint probabilities from the joint model
dataSim$pred.joint <- predict(mod, dataSim, type = "Biv")

#### predict marginals probabilities from the univariate models
dataSim$pred.naive <- predict(mod, dataSim, type = "naive")



######################################
#### Model estimation: Example 2  ####
######################################

%tau.eq1 <- seq(-0.45, 0.001, 0.05)
%tau.eq2 <- seq(-0.45, 0.001, 0.05)
%tau.eq1 <- ifelse(tau.eq1==0,0.001,tau.eq1)
%tau.eq2 <- ifelse(tau.eq2==0,0.001,tau.eq2)

%tau     <- expand.grid(tau.eq1 = tau.eq1, tau.eq2 = tau.eq2)
%tau$N.combinations <- 1:nrow(tau)



%out.N   <- list()


%# Gaussin Copula
%for (i in tau$N.combinations){
%  tryCatch({
%  tau.eq1 <-  tau[i,1]
% tau.eq2 <-  tau[i,2]
%    out.N[[i]]<- BivGEV(list(eq1,eq2), data = dataSim, BivD = "N", Model = "BivGEV", tau.eq1 = tau.eq1, tau.eq2 = tau.eq2)
%  }, error=function(e){paste("Tau combinations not supported")}
%  ,silent=FALSE)
%}



  }


\keyword{Bivariate Generalised Extreme Value Model}
