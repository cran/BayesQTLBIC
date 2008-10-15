\name{recalc.bicprobs}
\alias{recalc.bicprobs}
\title{Recalculation of posterior probabilities for bicreg.qtl  objects}
\description{Posterior probabilities for models, marginal probabilities
  for model sizes, and marginal probabilities for markers are
  re-calculated for new values of \code{delta} and/or prior
  probabilities}
\usage{recalc.bicprobs(obj, old.delta = attr(obj, "delta"), delta = 1, n = obj$n, 
    old.prior = attr(obj, "prior"), prior = 0.5, p.sg = 1) 
}
\arguments{
  \item{obj}{object of class bicreg.qtl}
  \item{old.delta}{value of \code{delta} previously used to calculate
    \code{obj}}
  \item{delta}{value of \code{delta} for which probabilities are
    desired}
  \item{n}{sample size used in the BIC calculation (equivalent to sample
    size before multiple imputations)}
  \item{old.prior}{(scalar) value of prior probability per marker used to
    calculate \code{obj}}
  \item{prior}{(scalar) value of prior probability per marker for which
    probabilities are desired}
  \item{p.sg}{proportion sampled in tails if selectively genotyping the population}
  }
  \details{Summary information in \code{obj} including BIC values,
    posterior probabilities for models, marginal probabilities for model
    size, and marginal probabilities for markers is recalculated with
    the specified new values of \code{delta} and \code{prior}
  }
\value{A list with components:
  \item{delta}{ new value of \code{delta} used}
  \item{prior}{ new value of \code{prior} used}
  \item{bic}{ new value of \code{bic} for each model},
  \item{postprob}{ new value of posterior probability for each model}
  \item{probne0}{ marginal probabilities for markers to be selected}
  \item{postprob.size}{ marginal probabilities for model size}
}
\author{R.D. Ball \email{rod.ball@AT@scionresearch.com}
}
\seealso{\code{\link{bicreg.qtl}}, \code{\link{summary.bicreg.qtl}}}
\examples{
  \dontrun{recalc.bicprobs(ex1n1200c1.bicreg,old.prior=0.2,prior=0.1)}
}
\keyword{regression}
\keyword{models}


  