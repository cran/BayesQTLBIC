\name{sim.bc.progeny}
\alias{sim.bc.progeny}
\title{Simulation of back-cross QTL progeny}
\description{Simulate a back-cross QTL progeny with given QTL
  and markers.}
\usage{sim.bc.progeny(n,Vp=NULL,map.pos,qtl.pos)}
\arguments{
  \item{n}{Number of progeny to simulate.}
  \item{Vp}{Vector of variances of QTL effects, as a proportion of total
    variance. Should sum to less than 1.}
  \item{map.pos}{Vector of map positions for markers, or a list
    with elements \code{map.pos} ( vector of map positions for markers)
    and \code{chrom} (vector of corresponding chromosome numbers).}
  \item{qtl.pos}{Vector of map positions for QTL, or a list
    with elements \code{qtl.pos} (vector of QTL positions)
    and \code{chrom} (vector of corresponding chromosome numbers).
    In case \code{Vp} is missing, \code{qtl.pos} should also contain elements
    \code{h2qs} (vector of QTL heritabilities) and
    \code{qtl.effect.signs} (vector of effect signs).}

}
\details{
Recombinations are simulated randomly along the genome for each progeny,
and used to generate marker and QTL genotypes. Random error with
variance \code{1 - sum(Vp)} is added to the QTL effects to give simulated trait values
\code{y}.
}
\value{
  \code{sim.bc.progeny} returns a list with elements \code{x} (matrix of marker
  values) and \code{y} (vector of trait values).
}
\author{R.D. Ball \email{rod.ball@scionresearch.com}}
\examples{
set.seed(1234)
ex1.map.pos <- seq(5,105,by=10)
qtldata200 <- sim.bc.progeny(n=200,Vp=c(0.1,0.2,0.3,0.15,0.25)/2,
      map.pos=list(chrom=rep(1:2,rep(length(ex1.map.pos),2)),
      pos=rep(ex1.map.pos,2)),qtl.pos=list(chrom=rep(1:2,c(3,2)),
                                pos=c(40,50,80,30,55)))
}
\keyword{regression}
\keyword{models}
