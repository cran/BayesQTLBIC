\name{impute.marker}
\alias{impute.marker}
\alias{impute.marker.matrix}
\title{Imputation of missing marker values in QTL progeny}
\description{Missing marker values are randomly imputed according to
  their probability, conditional on non-missing marker values.}

\usage{impute.marker(m, d, map.function="haldane", marker.values=1:2)

impute.marker.matrix(m, d, num.imputations = 1, map.function = "haldane", 
    marker.values = 1:2) }
\arguments{
  \item{m}{For \code{impute.marker}, a vector of marker haplotype values
    for a set of linked markers, 1 entry per marker in genome order.
    For \code{impute.marker.matrix}, a matrix with such vectors as rows
    corresponding to individual genotypes.}
  \item{d}{Inter-marker distances in cM; \code{d[i]} is the map distance
    from marker \code{i} to marker \code{i+1} 
%    (\eqn{i=1,\ldots,\mbox{length}(m)-1}{i=1,...,length(m)-1}).
     (\code{i=1,...,length(m)-1}).
  }
  \item{num.imputations}{Number of imputations to do.}
  \item{map.function}{Mapping function (Haldane or Kosambi) to use to
    estimate recombination probabilities.}
  \item{marker.values}{Vector of valid values for markers. Assuming
    diploid individuals and fully informative markers, the markers
    inherited from a given parent these would be one of two possible
    marker values.}
}

\details{For each imputation of each progeny, Missing marker values are
  randomly generated from their probability distribution, conditional on
  flanking markers where these exist, as follows. First, marker values
  not matching an entry of \code{marker.values} are replaced by missing
  values.  If there are missing values, one is randomly selected and
  sampled conditionally on any non-missing flanking marker values. The
  process is repeated (conditional on non-missing values, including
  originally missing values already sampled) until there are no
  missing values remaining. 
}

\value{\code{impute.marker} returns a vector with missing marker values
  randomly sampled conditional on non-missing values.  \code{impute.marker.matrix} returns 
  a matrix with missing marker values randomly sampled conditional on
  non-missing values in each row of \code{m}. If \code{num.imputations}
  is an integer greater than 1, \code{impute.marker.matrix} returns a matrix
  consisting of the \code{num.imputations} imputed matrices joined by row.}
\references{
 Ball, R. D. 2001: Bayesian methods for QTL mapping based on model selection: approximate analysis using the Bayesian Information
  Criterion. Genetics 159: 1351--1364.
}
\author{R.D. Ball \email{rod.ball@AT@scionresearch.com}}
\seealso{\code{\link{bicreg.qtl}}
}
\examples{
m <- c(1,1,2,2,1,NA,1,2,NA,1,1,NA)
d <- c(10,11,9,10,21,28,13,6,5,12,7)
impute.marker(m,d,map.function="haldane",marker.values=1:2)

# simulated backcross progeny with ~5\% missing marker genotypes
set.seed(1234)
marker.pos <- seq(5,105,by=10)
qtldata200 <- sim.bc.progeny(n=200,Vp=c(0.1,0.2,0.3,0.15,0.25)/2,
      map.pos=list(chrom=rep(1:2,rep(length(marker.pos),2)),
      pos=rep(marker.pos,2)),qtl.pos=list(chrom=rep(1:2,c(3,2)),
                                pos=c(40,50,80,30,55)))
x.mv <- qtldata200$x
num.missing <- round(0.05*length(qtldata200$x))
mv.pos <- sample(1:length(c(qtldata200$x)),size=num.missing,replace=FALSE)
x.mv[mv.pos] <- NA
c1.cols <- 1:11
x200c1.imp10 <- impute.marker.matrix(x.mv[,c1.cols],d=diff(marker.pos),
                                num.imputations=10)
c2.cols <- 12:22
x200c2.imp10 <- impute.marker.matrix(x.mv[,c2.cols],d=diff(marker.pos),
                                num.imputations=10)
x200 <- cbind(x200c1.imp10,x200c2.imp10)
}
\keyword{regression}
\keyword{models}
