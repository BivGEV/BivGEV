\name{summary.BivGEV}
\alias{summary.BivGEV}
%\alias{print.summary.BivGEV}

\title{
The function return some summaries from the \code{BivGEV} object.
}
\description{
The function use a fitted \code{BivGEV} object to return some summaries informations. %\code{summary.BivGEV} funcion is a wrapper from \code{summary.gjrm} in \code{GJRM}.
}
\usage{
\method{summary}{BivGEV}(object, n.sim = 1000, prob.lev = 0.05, ...)
}

\arguments{
  \item{object}{
A fitted \code{BivGEV} object.
}
 \item{n.sim}{
The simulated coefficient vectors from the posterior distribution of the estimated model parameters. The default number is 1,000.
}
  \item{prob.lev}{
Probability of the posterior distribution used for interval calculations.
}
\item{...}{Other arguments.}
}

\author{
Raffaella Calabrese - \email{raffaella.calabrese@ed.ac.uk};
Silvia Angela Osmetti - \email{silvia.osmetti@unicatt.it}; 
Mantainer: Luca Zanin - \email{luca.zanin@studio.unibo.it};
}

\seealso{
\code{\link{BivGEV}}
}


\keyword{Bivariate Generalised Extreme Value Model}

