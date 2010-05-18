\name{Unique.Rows}
\alias{Unique.Rows}
\title{Unique rows}
\description{Return indices of unique rows of a matrix}
\usage{
Unique.Rows(x)
}
\arguments{
  \item{x}{matrix}
}
\details{The \code{duplicated} function is applied to the rows of x}
\value{vector of indices of unique rows}
\author{R.D. Ball \email{rod.ball@scionresearch.com}}
\examples{
m <- matrix(c(1:3,1:3,2:4,5:7),nc=3,byrow=TRUE)
Unique.Rows(m)
m[Unique.Rows(m),]
}

