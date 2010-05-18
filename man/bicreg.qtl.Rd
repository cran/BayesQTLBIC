\name{bicreg.qtl}
\alias{bicreg.qtl}
\alias{bicreg2}
\title{Bayesian QTL analysis using the BIC criterion}
\description{Bayesian multi-locus QTL analysis based on Bayesian model
  selection in linear models using the BIC criterion to calculate
  approximate posterior probabilities for models.}
\usage{bicreg.qtl(x, y, wt = rep(1, length(y)), strict = TRUE, OR = 1000, 
    maxCol = 31, drop.factor.levels = TRUE, nvmax = 4, nbest = 10, 
    intercept = TRUE, do.occam = FALSE, n = length(y)/num.imputations, 
    num.imputations = 1, prior = 0.5, delta = 1, p.sg = 1, eval.markers = TRUE, 
    neval = NULL, keep.size = 1, method = c("regsubsets", "leaps")[1]) 
bicreg2(x, y, wt = rep(1, length(y)), strict = FALSE, OR = 1000, 
    maxCol = 31, drop.factor.levels = TRUE, nvmax = 4, nbest = 10, 
    intercept = TRUE, do.occam = TRUE, n = length(y)/num.imputations, 
    num.imputations = 1, prior = 0.5, delta = 1, p.sg = 1, eval.markers = FALSE, 
    neval = NULL, keep.size = -1, method = c("regsubsets", "leaps")[1]) 
}
\arguments{
  \item{x}{ a matrix of independent variables, based on marker
    genotypes, often from a single chromosome }
  \item{y}{ a vector of values for the dependent variable (trait values) }
  \item{wt}{ a vector of weights for regression }
  \item{strict}{ logical. \code{FALSE} returns all models whose posterior model probability is within a factor of \code{1/OR} of that of the best model. \code{TRUE} returns a more parsimonious set of models, where any model with a more likely submodel is eliminated.  }
  \item{OR}{ a number specifying the maximum ratio for excluding models in Occam's window }
  \item{maxCol}{ a number specifying the maximum number of columns in
    the design matrix (including the intercept) to be kept, i.e.
    maximum number of markers to include}
  \item{drop.factor.levels}{ logical. Indicates whether factor levels can be individually dropped in the stepwise procedure to reduce the number of columns in the design matrix, or if a factor can be dropped only in its entirety. }
  \item{nvmax}{ maximum number of variables in a model}
  \item{nbest}{ a value specifying the number of models of each size returned to bic.glm by the leaps algorithm.}
  \item{intercept}{ add an intercept term}
  \item{do.occam}{ do Occam's razor selection}
  \item{n}{ original sample size, before multiple imputations}
  \item{num.imputations}{number of imputations used to construct x, y}
  \item{prior}{ a vector or scalar specifying prior probabilities per
    marker for a QTL to be in the vicinity of the marker; generally
    proportional to the distance to flanking markers and total number of QTL
    expected genome. Defaults to 0.5 which is usually too high.}
  \item{delta}{ adjustment factor for the penalty term in the BIC criterion, 
    default is no adjustment (\code{delta=1}); (Cf. Broman and Speed 2002); not
    needed if using subjective prior probabilities and sample size is
    ample (\code{p.sg=1} and \code{n >= 100}; Ball 2007).
  }
  \item{p.sg}{ proportion genotyped (\code{p.sg/2} per tail), if selective
    genotyping is being used, default 1, corresponding to fully
    genotyped population}
  \item{eval.markers}{evaluate model averaged estimates for
    marker effects (effects of allelic substitution)}
  \item{neval}{ use \code{neval} top models on which to evaluate model averaged
    estimates of marker effects, default NULL, use all models}
  \item{keep.size}{ keep models up to this size regardless of Occam's
    razor criterion, e.g. to ensure the intercept only model is available for comparison}
  \item{method}{ choice of method, leaps or regsubsets}
}
\details{Provides posterior probabilities for linear models representing
  alternative QTL genetic architectures, which can be used for
  Bayesian inference of the  number of QTL and probabilities for QTL
  presence in a region. Provides Bayesian model averaged estimates for
  effects of QTL or effects of allelic substitution for markers which
  may be linked to QTL. Posterior probabilities are estimated from the
  BIC criterion combined with prior information, with adjustments for
  multiple imputation and selective genotyping. The posterior
  probability for model \eqn{M_i} is given by:

  \deqn{\mbox{Pr}(M_i) \propto \exp(-\mbox{BIC}_i/2) \times \pi(M_i)}{Pr(M_i) \\propto \exp(-BIC_i/2) \\times \pi(M_i)}

  where \eqn{\mbox{BIC}_i}{BIC_i} is the value of the BIC criterion and \eqn{\pi(M_i)} is 
  the prior probability for \eqn{M_i}.

  Missing marker values can be estimated by multiple imputation,
  conditional on flanking markers, using \code{impute.markers},  and the
  imputed data used as x ,y.

  For selectively genotyped populations (Darvasi and Soller 1992) an
  adjustment is made to the BIC criterion. Asymptotic convergence is good for fully genotyped
  families  with \eqn{n\geq 100}{n >= 100} progeny but requires larger sample sizes for
  smaller proportions of the tails (\code{p.sg})  genotyped.
}
\value{\code{bicreg.qtl} returns an object of class \code{bicreg.qtl}

The function \code{summary} is used to print a summary of the results. 

An object of class \code{bicreg.qtl} inherits from class \code{bicreg}
and is a list containing at least the following components/attributes:

  \item{postprob}{the posterior probabilities of the models selected}
  \item{namesx}{the names of the variables}
  \item{label}{labels identifying the models selected}
  \item{r2}{R2 values for the models}
  \item{bic}{values of BIC for the models}
  \item{size}{the number of independent variables in each of the models}
  \item{which}{a logical matrix with one row per model and one column per variable indicating whether that variable is in the model}
  \item{probne0}{the posterior probability that each variable is non-zero (in percent)}
\item{n}{the sample size before multiple imputation}
\item{postprob.size}{the marginal posterior probabilities for model sizes}
  \item{postmean}{the posterior mean of each coefficient (from model averaging)}
  \item{postsd}{the posterior standard deviation of each coefficient (from model averaging) }
  \item{condpostmean}{the posterior mean of each coefficient conditional on the variable being included in the model}
  \item{condpostsd}{the posterior standard deviation of each coefficient conditional on the variable being included in the model}
  \item{ols}{matrix with one row per model and one column per variable giving the OLS estimate of each coefficient for each model}
  \item{se}{matrix with one row per model and one column per variable giving the standard error of each coefficient for each model}
  \item{reduced}{a logical indicating whether any variables were dropped before model averaging}
  \item{dropped}{a vector containing the names of those variables dropped before model averaging}
  \item{call}{the matched call that created the bicreg object}
  and containts the following attributes:
  \item{intercept}{if an intercept term was added}
  \item{num.imputations}{the number of multiple imputations assumed}
  \item{p.sg}{the proportion genotyped (\code{p.sg/2} per tail)}
  \item{delta}{the value of delta used}
  }
\references{
  Ball, R. D. 2001: Bayesian methods for QTL mapping based on model selection: approximate analysis using the Bayesian Information
  Criterion. Genetics 159: 1351--1364.

  Ball, R.D. 2007:   Quantifying evidence for candidate gene
  polymorphisms---Bayesian analysis combining sequence-specific and QTL
  co-location information. Genetics 177: 2399--2416. 

  DeSilva, H.N., and Ball, R.D. 2007:  Linkage disequilibrium mapping concepts.
  Chapter 7, pp103--132 In: Association Mapping in Plants, N.C. Oraguzie,
  E.H.A. Rikkerink, S.E. Gardiner, and H.N. DeSilva (Editors), Springer,
  New York.

  Ball, R.D. 2007:  Statistical analysis and experimental design.
  Chapter 8, pp133--196 In: Association Mapping in Plants, N.C. Oraguzie,
  E.H.A. Rikkerink, S.E. Gardiner, and H.N. DeSilva (Editors), Springer,
  New York.

  Bogdan, M., Ghosh, J. K.,  and Doerge, R. W. 2004: Modifying the Schwarz
  Bayesian information criterion to locate multiple interacting
  quantitative trait loci. Genetics 167: 989--999.

  Broman, K.W. and Speed, T.P. 2002: A model selection approach for the
  identification of quantitative trait loci in experimental crosses
  (with discussion). J Roy Stat Soc B 64: 641--656, 731--775.

  Darvasi, A. and Soller, M. 1992: Selective genotyping for determination
  of linkage between a locus and a quantitative trait locus.
  Theoretical and Applied Genetics 85: 353--359.

  Raftery, A. E. 1995:  Bayesian model selection in social research (with 
  Discussion). Sociological Methodology 1995 (Peter V. Marsden, ed.),
  pp. 111-196, Cambridge, Mass.: Blackwells.
  
}
\author{R.D. Ball, 
%Scion (New Zealand Forest Research Institute Limited);
  (\email{rod.ball@scionresearch.com}),
  based on \code{bicreg} from the \code{BMA} package by Adrian Raftery,
  Chris Volinski, and Ian Painter. }
\seealso{\code{bicreg}, 
  \code{\link{summary.bicreg.qtl}},
  \code{\link{DS.gamma}}, \code{\link{recalc.bicprobs}}
  \code{\link{impute.marker}} }
\examples{
# simulated backcross progeny 
set.seed(1234)
ex1.marker.pos <- seq(5,105,by=10)
chrom <- rep(1:2,rep(length(ex1.marker.pos),2))
ex1.qtldata <- sim.bc.progeny(n=1200,Vp=c(0.1,0.2,0.3,0.15,0.25)/2,
           map.pos=list(chrom=rep(1:2,rep(length(ex1.marker.pos),2)),
           pos=rep(ex1.marker.pos,2)),qtl.pos=list(chrom=rep(1:2,c(3,2)),
                                     pos=c(40,50,80,30,55)))
ex1n1200c1.bicreg <- bicreg.qtl(x=ex1.qtldata$x[,chrom==1],y=ex1.qtldata$y,OR=1000,
                                nbest=10,nvmax=5,prior=0.2,keep.size=1)
# 23 models account for 99\% of the probability
cumsum(ex1n1200c1.bicreg$postprob)
# 2 QTL in coupling at 40,50cM can't be resolved
summary(ex1n1200c1.bicreg,nbest=23)
}
\keyword{regression}
\keyword{models}
