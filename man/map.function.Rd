\name{map.function.rec}
\alias{map.function.rec}
\alias{map.function.dist}
\title{Map functions}
\description{Map functions to convert between recombination and physical
  map distances for the Haldane and Kosambi mapping functions}
\usage{
map.function.rec(d,method="haldane")
map.function.dist(r,method="haldane")
}
\arguments{
  \item{d}{map distance in Morgans}
  \item{r}{recombination distance}
  \item{method}{choice of mapping function}
}
\details{The conversion between physical distance \code{d} (Morgans)
  and recombination probabilities is given by:
  \deqn{r = \frac{1}{2} (1 - \exp(-2d))}{r = 1/2*(1-exp(-2d))}
  \deqn{d = -\frac{1}{2} \log(1-2r)}{d = -1/2*log(1-2r)}
for the Haldane mapping function, and by:
  \deqn{r = \frac{1}{2}\frac{1- \exp(-4d)}{1+\exp(-4d)}}{r = 1/2*(1- \exp(-4d))/(1+\exp(-4d))}
  \deqn{d =\frac{1}{4}\log\left(\frac{1+2r}{1-2r}\right)}{d = 1/4*log((1+2*r)/(1-2*r))}
for the Kosambi mapping function.
}
\value{the value of \code{r} or \code{d}}
\references{
Liu, B. 1998: Statistical Genomics: Linkage, Mapping and QTL Analysis.
Boca Raton, FL: CRC Press LLC.
}
\author{R.D. Ball \email{rod.ball@AT@scionresearch.com}}
\examples{
 map.function.rec(d=c(0,0.01,0.1,0.2,1,5),method="haldane")
 map.function.rec(d=c(0.01,0.1,0.2,1,5),method="kosambi")
 map.function.dist(r=c(0,0.01,0.1,0.2,0.49,0.499,0.5),method="haldane")
 map.function.dist(r=c(0,0.01,0.1,0.2,0.49,0.499,0.5),method="kosambi")
}
\keyword{models}
