\name{sample.bicreg.qtl.models}
\alias{sample.bicreg.qtl.models}
\title{Sampling models for BIC analyses of multiple chromosomes}
\description{Models are sampled from each of a set of runs with
  different sets of independent variables (typically
  corresponding to multiple chromosomes), according to their posterior
  probabilities.}
\usage{sample.bicreg.qtl.models(chrom.fits,nsim,maxtries=10)}
\arguments{
  \item{chrom.fits}{list of \code{bicreg.qtl} objects from separate
    fits by chromosome or genomic region}
  \item{nsim}{number of models to sample}
  \item{maxtries}{maximum number of retries to get nsim unique models}
}
\details{Each of \code{nsim} combined models is obtained by randomly sampling one model from each
  chromosome according to its posterior probability, and combining the
  \code{x}-variables from the sampled models.
}
\value{A list of models represented as a matrix (similar to the
  \code{which} matrix returned by \code{\link{bicreg.qtl}},
   whose \code{(i,j)} element is \code{TRUE} if the \code{i}th
  sampled model contains the \code{j}th variable}
\references{
  Ball, R. D. 2001: Bayesian methods for QTL mapping based on model selection: approximate analysis using the Bayesian Information
  Criterion. Genetics 159: 1351--1364.}
\author{R.D. Ball, (\email{rod.ball@AT@scionresearch.com})}
\seealso{\code{\link{bicreg.qtl}},\code{\link{bicreg.models}}}
\examples{
  \dontrun{mWhich200 <- sample.bicreg.qtl.models(chrom.fits,nsim=200)}}
\keyword{regression}
\keyword{models}
