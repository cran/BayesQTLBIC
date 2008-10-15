\name{bicreg.models}
\alias{bicreg.models}
\title{BIC analysis for a fixed set of models}
\description{\code{bicreg.models} evaluates posterior probabilities based on
  the BIC criterion and prior probabilities for a fixed set of models,
  typically obtained by sampling from sets of models from separate analyses 
  of individual chromosomes.
}
\usage{bicreg.models(x,y,wt = rep(1, length(y)),which,intercept=TRUE,add.null.model=TRUE,
                          n=length(y)/num.imputations,num.imputations=1,delta=1,
                          p.sg=1,prior=0.5,eval.markers=TRUE,neval=NULL)
}
\arguments{
  \item{x}{ Matrix of independent variables, based on marker
    genotypes, often from a single chromosome. }
  \item{y}{ Vector of values for the dependent variable (trait values). }
  \item{wt}{ Vector of weights for regression. }
  \item{which}{Matrix of logical values corresponding to a set of
    models, the (i,j) element is TRUE if and only if the jth variable is
    selected in the ith model. }
  \item{intercept}{ Add an intercept term.}
  \item{add.null.model}{ Add the NULL model.}
  \item{n}{ Original sample size, before multiple imputations.}
  \item{num.imputations}{Number of imputations used to construct x, y.}
  \item{prior}{ Vector or scalar specifying prior probabilities per
    marker for a QTL to be in the vicinity of the marker; generally
    proportional to the distance to flanking markers and total number of QTL
    expected genome. Defaults to 0.5 which is usually too high.}
  \item{delta}{ Adjustment factor for the penalty term in the BIC criterion, 
    default is no adjustment delta=1; (Cf. Broman and Speed 2002); not
    needed if using subjective prior probabilities and sample size is
    ample (\code{p.sg=1} and \eqn{n \ge 100}; Ball 2007).
  }
  \item{p.sg}{ Proportion \code{p.sg/2} of each tail is genotyped if selective
    genotyping is being used; default 1, corresponding to fully
    genotyped population.}
  \item{eval.markers}{ Evaluate model averaged estimates for
    marker effects (effects of allelic substitution).}
  \item{neval}{ Use \code{neval} top models on which to evaluate model averaged
    estimates of marker effects, default NULL, use all models.}
}
\details{Provides posterior probabilities for a fixed set of linear
  models representing alternative QTL genetic architectures. Provides
  Bayesian model averaged estimates for effects of QTL or effects of allelic
  substitution for markers which may be linked to QTL.}
\value{\code{bicreg.models}  returns an object of class
  \code{bicreg.qtl}.}
\references{
  Ball, R. D. 2001: Bayesian methods for QTL mapping based on model selection: approximate analysis using the Bayesian Information
  Criterion. Genetics 159: 1351--1364.}
\author{R.D. Ball, (\email{rod.ball@AT@scionresearch.com})}
\seealso{\code{\link{bicreg.qtl}},\code{\link{sample.bicreg.qtl.models}}}
\examples{
\dontrun{
data(ex3n300a.data)
chrom <- rep(1:12,rep(16,12))
marker <- rep(1:16,12)
x <- sapply(ex3n300a.data$Markers,c)
y <- ex3n300a.data$Trait$t1
nchrom <- length(sort(chrom.levels <- unique(chrom)))

quick.demo <- TRUE
if(quick.demo){
  nc <- 2; nsim <- 20;x <- x[,chrom \%in\% 1:2];
  chrom <- chrom[chrom \%in\% 1:2]
}else{nc <- 12; nsim <- 200}
chrom.fits <- list()
for(ii in seq(along=chrom.levels[1:nc])){
  cat(paste("*** chromosome",ii,"***","\n"))
  ci <- chrom.levels[ii]
  chrom.sel <- chrom==ci
  chrom.fits[[ii]] <- bicreg.qtl(x[,chrom.sel],y, prior=0.1,nbest=20,nvmax=3)
}
mWhich <- sample.bicreg.qtl.models(chrom.fits,nsim=nsim)
mres <- bicreg.models(x=x,y=y,which=mWhich,prior=0.1)
summary(mres,nbest=38,min.marker.prob=0.05)
}
}
\keyword{regression}
\keyword{models}
