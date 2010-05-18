\name{DS.gamma}
\alias{DS.gamma}
\title{Darvasi and Soller correction factor gamma}
\description{Calculate Darvasi and Soller correction factor gamma for
  QTL mapping in selectively genotyped populations}
\usage{DS.gamma(p)}
\arguments{
  \item{p}{proportion sampled (\eqn{p/2} per sampled tail)} 
}
\details{For a given number of individuals genotyped, power of detection of 
  QTL for a continuous trait (assumed normally distributed) is increased
  if individuals are sampled from the tails of the  phenotypic distribution. 
  For small gene effects, if proportion p/2 is sampled from each tail of
  the phenotypic distribution, naive estimates of marker effects are
  inflated by the factor \eqn{\gamma_p} (to first order in the size of
  gene effects) where:
  \deqn{\gamma_p = 1 +
    \Phi^{-1}(1-p/2)\frac{\phi(\Phi^{-1}(p/2))}{p/2}.}{
     gamma_p = 1 + qnorm((1 - p/2)*dnorm(qnorm(p/2))/(p/2))
    }
}
\value{the value of \eqn{\gamma_p}{gamma_p}}
\references{
  Darvasi, A. and Soller, M. 1992: Selective genotyping for determination
  of linkage between a locus and a quantitative trait locus.
  Theoretical and Applied Genetics 85: 353--359.
}
\author{Rod Ball (\email{rod.ball@scionresearch.com})}
\examples{
DS.gamma(p = c(0.05,0.1,0.25))
}
\keyword{models}