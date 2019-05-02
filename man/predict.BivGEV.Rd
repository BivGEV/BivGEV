\name{predict.BivGEV}
\alias{predict.BivGEV}

\title{
This function produces predictions for univariate and joint GEV models.
}
\description{
The function takes the fitted \code{BivGEV} object and produces univariate and joint predicted probabilities for a set of values (new or used for the training) of the model.
}
\usage{
\method{predict}{BivGEV}(object, newdata, type = "joint")
}

\arguments{
  \item{object}{
A fitted \code{BivGEV} object.
}
  \item{newdata}{
The data frame with the covariates of the model.
}
  \item{type}{
The possible values are "naive" (univariate estimates) and "joint" (bivariate estimates).
}
}
\details{
\itemize{
\item The function allows to compute the marginal probabilities of each equation in the assumption that the two models are estimated separatelly by selecting the option type = "naive". Specifically, the option return the marginal probabilities \eqn{Pr(Y1)}=1 and \eqn{Pr(Y2)}=1 named \emph{p1.naive} and \emph{p2.naive}, respectively. 

\item Selecting the option type = "joint" the function compute the probabilities from the joint model. Specifically, the option return the marginal probabilities \eqn{Pr(Y1)}=1 and \eqn{Pr(Y2)}=1 and the joint probabilities \eqn{Pr(Y1=1,Y2=1)}, \eqn{Pr(Y1=1,Y2=0)}, \eqn{Pr(Y1=0,Y2=1)}, \eqn{Pr(Y1=0,Y2=0)} named \emph{p1}, \emph{p2}, \emph{p11}, \emph{p10}, \emph{p01}, and \emph{p00}, respectively.
}
}

\author{
Raffaella Calabrese - \email{raffaella.calabrese@ed.ac.uk};
Silvia Angela Osmetti - \email{silvia.osmetti@unicatt.it}; 
Luca Zanin - \email{luca.zanin@studio.unibo.it};
}

\seealso{
\code{\link{BivGEV}}
}


\keyword{Bivariate Generalised Extreme Value Model}
