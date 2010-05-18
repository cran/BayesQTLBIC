#require(leaps)

bicreg.qtl <-function(x, y, wt = rep(1, length(y)), strict = TRUE, OR = 1000, maxCol=31,
         drop.factor.levels=TRUE, nvmax=4, nbest=10, intercept=TRUE, do.occam=FALSE,
         n = length(y)/num.imputations, num.imputations = 1, prior = 0.5, delta=1,
	 p.sg=1, eval.markers=TRUE, neval=NULL,keep.size=1,
         method=c("regsubsets","leaps")[1]){
  dropcols <- function(x, y, wt, maxCols = 31) {
    x1.ldf <- data.frame(x, y = y)
    temp.wt <- wt
    lm.out <- lm(y ~ ., data = x1.ldf, weights = temp.wt)
    form.vars <- all.vars(formula(lm.out))[-1]
    any.dropped <- FALSE
    dropped.which <- NULL
    while (length(lm.out$coefficients) > maxCol) {
      any.dropped <- TRUE
      droplm <- drop1(lm.out, test = "none")
      dropped <- row.names(droplm)[which.min(droplm$RSS[-1]) + 
                                   1]
      dropped.index <- match(dropped, form.vars)
      form.vars <- form.vars[-dropped.index]
      formla <- formula(paste("y", "~", paste(form.vars, 
                                              collapse = " + "), sep = " "))
      lm.out <- lm(formla, data = x1.ldf, weights = temp.wt)
      dropped.which <- c(dropped.which, dropped)
    }
    new.var.names <- names(lm.out$coefficients)
    return(list(mm = model.matrix(lm.out)[, -1, drop = FALSE], 
                any.dropped = any.dropped, dropped = dropped.which, 
                var.names = new.var.names))
  }
  cl <- match.call()
  x <- data.frame(x)
  if (is.null(dimnames(x))) 
    dimnames(x) <- list(NULL, paste("X", 1:ncol(x), sep = ""))
  y <- as.numeric(y)
  options(contrasts = c("contr.treatment", "contr.treatment"))
  xnames <- dimnames(x)[[2]]
  x2 <- na.omit(data.frame(x))
  used <- match(row.names(data.frame(x)), row.names(x2))
  omitted <- seq(nrow(x))[is.na(used)]
  if (length(omitted) > 0) {
    wt <- wt[-omitted]
    x <- x2
    y <- y[-omitted]
    warning(paste("There were ", length(omitted), "records deleted due to NA's"))
  }
  if (drop.factor.levels) {
    cdf <- cbind.data.frame(y = y, x)
    mm <- model.matrix(formula(cdf), data = cdf)[, -1, drop = FALSE]
    x <- mm
  }
  xx <- dropcols(x, y, wt, maxCol)
  xnames <- xx$var.names[-1]
  x <- xx$mm
  reduced <- xx$any.dropped
  dropped <- NULL
  if (reduced) 
    dropped <- xx$dropped
  nvar <- length(x[1, ])
  if (nvar > 2) {
    if(!method  %in% c("leaps","regsubsets"))stop("method should be either \"leaps\" or \"regsubsets\"")
    if(method=="leaps"){
      a <- leaps(x, y, wt = wt, int=intercept, method = "r2", names = dimnames(x)[[2]], 
                 strictly.compatible = FALSE, nbest = nbest)
      a$r2 <- pmin(pmax(0, a$r2), 0.999)
#      x.lm <- cbind.data.frame(y = y, as.data.frame(x[, a$which[2, 
#                                 , drop = FALSE]]), w = wt)
#      lm.fix <- lm(y ~ . -w, weights = w, data = x.lm)
# RDB R CMD check doesn't find w
# Splus problem required w in the dataframe ?
      x.lm <- cbind.data.frame(y = y, as.data.frame(x[, a$which[2, 
                                 , drop = FALSE]]))
      lm.fix <- lm(y ~ . , weights = wt, data = x.lm)
      r2.fix <- summary(lm.fix)$r.sq
      N <- ncol(x)
      magic <- N * log(1 - a$r2[2]) - N * log(1 - r2.fix)
      a$r2 <- 1 - (1 - a$r2) * exp(-magic/N)
      r2 <- round(c(0, a$r2) * 100, 3)
      size <- c(1, a$size)
      which <- rbind(rep(FALSE, ncol(x)), a$which)
      templabs <- t(matrix(rep(colnames(which), times = nrow(which)), 
                           ncol = nrow(which)))
      templabs[!which] <- ""
      label <- apply(templabs, 1, paste, collapse = "")
      label[1] <- "NULL"
    }else{ # method=="regsubsets
       search.method <- "exhaustive"
       a1 <- regsubsets(x=x,y=y,weights=wt,nbest=nbest,nvmax=nvmax,intercept=intercept, method=search.method,really.big=T)
       s1 <- summary(a1,matrix=TRUE,matrix.logical=TRUE)
       a <- list(which=s1$which)
       a$size <- 1 + c(0,as.numeric(rownames(a$which)))
       a$which <- rbind(c(TRUE,rep(FALSE,ncol(x))),a$which)[,-1]
       rownames(a$which)[1] <- "0"
       a$r2 <- s1$rsq
       a$label <- colnames(a$which)
       r2 <- round(c(0, a$r2) * 100, 3)
       size <- a$size
       label <- a$label
       which <- a$which
      }
  }else{ # 2 or less variables
    r2 <- bic <- NULL
    nmod <- switch(ncol(x), 2, 4)
    bic <- label <- rep(0, nmod)
    model.fits <- as.list(rep(0, nmod))
    which <- matrix(c(FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, 
                      TRUE, TRUE), nmod, nmod/2)
    size <- c(1, 2, 2, 3)[1:nmod]
    sep <- if (all(nchar(dimnames(x)[[2]]) == 1)) 
      ""
    else ","
    for (k in 1:nmod) {
      if (k == 1) {
        label[k] <- "NULL"
        lm1 <- lm(y ~ 1, w = wt)
      }
      else {
        label[k] <- paste(dimnames(x)[[2]][which[k, ]], 
                          collapse = sep)
        x.lm <- cbind.data.frame(y = y, x = x[, which[k, 
                                          , drop = FALSE]], wt = wt)
        lm1 <- lm(y ~ . - wt, data = x.lm, w = wt)
      }
      r2[k] <- summary(lm1)$r.sq * 100
    }
  }

  ## RB: adjust for multiple imputations
#browser()
  n <- length(y)/num.imputations
  bic <- n * log(1 - r2/100) + (size - 1) * log(n)

  ## RB: incorporate adjustment factor delta
  if(!missing(delta) && delta != 1){
    bic <- bic.delta <- n * log(1 - r2/100) + delta*(size - 1) * log(n)
  }
  ## RB: adjustment for selective genotyping
  if(!missing(p.sg) && p.sg < 1 && p.sg > 0){
    gammap <- DS.gamma(p.sg)
    log.gammap <- log(gammap)
    bic <- bic + log.gammap*(size - 1)
  }
  ## RB: adjust for prior
  nmod <- length(bic)
  if(!missing(prior)) {
    if(length(prior) == 1)
      prior <- rep(prior, nvar)
    log.prior <- log(prior)
    logc.prior <- log(1 - prior)
    log.prior.probs <- rep(NA, nmod)
    for(ii in 1:nmod) {
      log.prior.probs[ii] <- sum(which[ii,  ] * log.prior + (1 - which[ii,  ]) * logc.prior)
    }
    bic <- bic - 2 * log.prior.probs
  }
  if(do.occam){
    occam <- bic - min(bic) < 2 * log(OR)
    if(!is.null(keep.size)) occam | size <= keep.size
  }else{
    occam <- rep(TRUE,length(bic))
  }
  r2 <- r2[occam]
  size <- size[occam]
  label <- label[occam]
  which <- which[occam, , drop = FALSE]
  bic <- bic[occam]
  mbic <- bic - max(bic)
  postprob <- exp(-0.5 * mbic)/sum(exp(-0.5 * mbic))
  postprob[is.na(postprob)] <- 1
  order.bic <- order(bic, size)
  r2 <- r2[order.bic]
  size <- size[order.bic]
  label <- label[order.bic]
  which <- which[order.bic, , drop = FALSE]
  bic <- bic[order.bic]
  postprob <- postprob[order.bic]
  if (strict & do.occam) {
    nmod <- length(bic)
    if (nmod > 1) {
      occam <- rep(TRUE, nmod)
      for (k in (2:nmod)) {
        for (j in (1:(k - 1))) {
          which.diff <- which[k, ] - which[j, ]
          if (all(which.diff >= 0)) 
            occam[k] <- FALSE
        }
      }
      r2 <- r2[occam]
      size <- size[occam]
      label <- label[occam]
      nmod <- sum(occam)
      which <- which[occam, , drop = FALSE]
      bic <- bic[occam]
      postprob <- postprob[occam]
      postprob <- postprob/sum(postprob)
    }
  }
  probne0 <- round(100 * t(which) %*% as.matrix(postprob), 1)
  nmod <- length(bic)
  model.fits <- as.list(rep(0, nmod))
  for (k in (1:nmod)) {
    if (sum(which[k, ]) != 0) {
      model.fits[[k]] <- ls.print(lsfit(x[, which[k, ], 
                                          drop = FALSE], y, wt = wt), print.it = FALSE)$coef.table[[1]]
    }
    else model.fits[[k]] <- ls.print(lsfit(rep(1, length(y)), 
                                           y, wt = wt, int = FALSE), print.it = FALSE)$coef.table[[1]]
  }
  Ebi <- rep(0, (nvar + 1))
  SDbi <- rep(0, (nvar + 1))
  CEbi <- Ebi
  CSDbi <- SDbi
  EbiMk <- matrix(rep(0, nmod * (nvar + 1)), nrow = nmod)
  sebiMk <- matrix(rep(0, nmod * (nvar + 1)), nrow = nmod)
  for (i in 1:(nvar + 1)) {
    if ((i == 1) || (sum(which[, (i - 1)] != 0))) {
      for (k in (1:nmod)) {
        if ((i == 1) || (which[k, (i - 1)] == TRUE)) {
          if (i == 1) 
            pos <- 1
          else pos <- 1 + sum(which[k, (1:(i - 1))])
          EbiMk[k, i] <- model.fits[[k]][pos, 1]
          sebiMk[k, i] <- model.fits[[k]][pos, 2]
        }
      }
      Ebi[i] <- as.numeric(sum(postprob * EbiMk[, i]))
      SDbi[i] <- sqrt(postprob %*% (sebiMk[, i]^2) + postprob %*% 
                      ((EbiMk[, i] - Ebi[i])^2))
      if (i == 1) {
        CEbi[i] <- Ebi[i]
        CSDbi[i] <- SDbi[i]
      }
      else {
        sel <- which[, i - 1]
        cpp <- postprob[sel]/sum(postprob[sel])
        CEbi[i] <- as.numeric(sum(cpp * EbiMk[sel, i]))
        CSDbi[i] <- sqrt(cpp %*% (sebiMk[sel, i]^2) + 
                         cpp %*% ((EbiMk[sel, i] - CEbi[i])^2))
      }
    }
  }
  dimnames(which) <- list(NULL, dimnames(x)[[2]])
  dimnames(EbiMk) <- dimnames(sebiMk) <- list(NULL, c("Int", 
                                                      dimnames(x)[[2]]))
  ## RB: evaluate effects of each marker (models of size 1) against other models
  ## RB: effects of allelic substitution
  if(eval.markers){
    if(is.null(neval) || neval==0)neval <- nmod
    d <- matrix(NA,nr=neval,nc=ncol(x))
    for(ii in 1:neval){
      if(size[ii]>1){
        mui <- model.fits[[ii]][1,1]
        ypi <- as.matrix(x)[,which[ii,],drop=F] %*% model.fits[[ii]][-1,1]
        for(jj in 1:ncol(x)){
          d[ii,jj] <- lm(ypi ~ x[,jj])$coef[2]
        }
      }else d[ii,] <- rep(0,ncol(d))
    }
    d.av <- c(postprob %*% d)
  }
  ## RB: calculate marginal probabilities for model size
  ## done in other BMA functions?
  max.size <- max(size)
  mprob.size <- tapply(c(postprob,rep(0,max.size)),c(size,1:max.size)-1,
                       sum)
  if(is.matrix(probne0)){
    marker.names <- dimnames(probne0)[[1]]
    probne0 <- c(probne0)
    names(probne0) <- marker.names
  }

  result <- list(postprob = postprob, namesx = xnames, label = label, 
                 r2 = r2, bic = bic, size = (size - 1), which = which, 
                 probne0 = c(probne0),
                 n = n, ## RB added for recalc.bicprobs
                 postprob.size=mprob.size, ## RB added for compatibility with our functions
                 postmean = Ebi, postsd = SDbi, 
                 condpostmean = CEbi, condpostsd = CSDbi, ols = EbiMk, 
                 mle = EbiMk, se = sebiMk, reduced = reduced, dropped = dropped, 
                 call = cl, n.models = length(postprob), n.vars = length(probne0))
  if(eval.markers)result$d.av <- d.av
  if(length(prior)==1 || sqrt(var(prior)) < 1e-6*mean(prior))attr(result,"prior") <- prior[1] else attr(result,"prior") <- prior
  attr(result,"intercept") <- intercept
  attr(result,"num.imputations") <- num.imputations
  attr(result,"p.sg") <- p.sg
  attr(result,"delta") <- delta
  class(result) <- c("bicreg.qtl","bicreg")
  result
}

bicreg2 <- bicreg.qtl
formals(bicreg2)$do.occam <- TRUE
formals(bicreg2)$strict <- FALSE
formals(bicreg2)$eval.markers <- FALSE
formals(bicreg2)$keep.size <- -1

args(bicreg.qtl)
args(bicreg2)

sample.bicreg.qtl.models <-
  function (chrom.fits, nsim, maxtries=10) 
# chrom.fits, a list of fits from bicreg.qtl
# sample nsim random model with independently for each chromosome, 
# with replacement according to posterior probabilities
# non-duplicated rows returned
{
  res <- sample.bicreg.qtl.models1(chrom.fits,nsim)
  res <- res[Unique.Rows(res),]
  nr <- nrow(res)
  ntries <- 1
  while(nr < nsim && ntries < maxtries){ 
    res2 <- sample.bicreg.qtl.models1(chrom.fits,nsim)
    res <- rbind(res,res2)
    res <- res[Unique.Rows(res),]
    nr <- nrow(res)
  }
  res[1:nsim,]
}


sample.bicreg.qtl.models1 <- function(chrom.fits,nsim){
# chrom.fits, a list of fits from bicreg.qtl
# sample nsim random model with independently for each chromosome, 
# with replacement according to posterior probabilities
# rows may be duplicated
  l <- lapply(chrom.fits, function(fitii){
#browser()
   pprobii <- fitii$postprob
   nmodii <- length(pprobii)
   fitii$which[sample(1:nmodii,size=nsim,prob=pprobii,replace=TRUE),]})
 do.call(cbind,l)
}

sample.bicreg.qtl.models.old <- function(chrom.fits,nsim){
# chrom.fits, a list of fits from bicreg.qtl
# sample nsim random model with independently for each chromosome, 
# with replacement according to posterior probabilities

  l <- lapply(chrom.fits, function(fitii){
#browser()
   pprobii <- fitii$postprob
   nmodii <- length(pprobii)
   fitii$which[sample(1:nmodii,size=nsim,prob=pprobii,replace=TRUE),]})
if(F){
nchrom <- length(chrom.fits)
#browser()
#length(u)/(nchrom*nsim)
u <- sapply(l, function(m)c(t(m)))
u1 <- array(c(u),dim=c(length(u)/(nchrom*nsim),nsim,nchrom))
mWhich <- t(matrix(aperm(u1,c(1,3,2)),ncol=nsim))
mWhich
}

do.call(cbind,l)
}


# same as bicreg but with a fixed set of models

bicreg.models <-
  function (x, y, wt = rep(1, length(y)), which, intercept = TRUE, 
            add.null.model = TRUE, n = length(y)/num.imputations, num.imputations = 1, 
            delta = 1, p.sg = 1, prior = 0.5, eval.markers = TRUE, neval = NULL) 
{
  nvar <- length(x[1, ])
  if (!(is.matrix(which) & nrow(which) > 0)) 
    stop("need at least 1 row in which matrix")
  nmod <- nrow(which)
  model.fits <- as.list(rep(0, nmod))
  a <- list()
  a$size <- a$r2 <- numeric(nrow(which))
  a$which <- which
  for (ii in 1:nrow(which)) {
    if (sum(which[ii, ]) != 0) {
      summii <- ls.print(lsfit(x[, which[ii, ], drop = FALSE], 
                               y, wt = wt), print.it = FALSE)
      model.fits[[ii]] <- summii$coef.table[[1]]
      a$r2[ii] <- as.numeric(summii$summary[1, "R Squared"])

    }
    else {
#       gives wrong R^2 for null model
#      summii <- ls.print(lsfit(rep(1, length(y)), y, wt = wt, 
#                               int = FALSE), print.it = FALSE)
      summii <- summary(lm(y~1, weight=wt))
      model.fits[[ii]] <- coef(summii)
      a$r2[ii] <- summii$r.squared
    }
    a$size[ii] <- sum(which[ii, ])
  }
  #browser()
  if (add.null.model && min(a$size) > 0) {
    size <- 1 + c(0, a$size)
    r2 <- round(c(0, a$r2) * 100, 3)
    which <- rbind(rep(FALSE, ncol(a$which)), a$which)
    model.fits <- c(list(ls.print(lsfit(rep(1, length(y)), 
                                        y, wt = wt, int = FALSE), print.it = FALSE)$coef.table[[1]]), 
                    model.fits)
  } else {
    size <- 1 + a$size
    r2 <- round(a$r2 * 100, 3)
  }
  n <- length(y)/num.imputations
  bic <- n * log(1 - r2/100) + (size - 1) * log(n)
  if (!missing(delta) && delta != 1) {
    bic <- bic.delta <- n * log(1 - r2/100) + delta * (size - 
                                                       1) * log(n)
  }
  if (!missing(p.sg) && p.sg < 1 && p.sg > 0) {
    gammap <- DS.gamma(p.sg)
    log.gammap <- log(gammap)
    bic <- bic + 2 * log.gammap * (size - 1)
  }
  nmod <- length(bic)
  if (!missing(prior)) {
    if (length(prior) == 1) 
      prior <- rep(prior, nvar)
    log.prior <- log(prior)
    logc.prior <- log(1 - prior)
    log.prior.probs <- rep(NA, nmod)
    for (ii in 1:nmod) {
      log.prior.probs[ii] <- sum(which[ii, ] * log.prior + 
                                 (1 - which[ii, ]) * logc.prior)
    }
    bic <- bic - 2 * log.prior.probs
  }
  mbic <- bic - max(bic)
  postprob <- exp(-0.5 * mbic)/sum(exp(-0.5 * mbic))
  postprob[is.na(postprob)] <- 1
  order.bic <- order(bic, size)
  r2 <- r2[order.bic]
  size <- size[order.bic]
  which <- which[order.bic, , drop = FALSE]
  model.fits <- model.fits[order.bic]
  bic <- bic[order.bic]
  postprob <- postprob[order.bic]
  probne0 <- round(100 * t(which) %*% as.matrix(postprob), 1)
  nmod <- length(bic)
  cat("before calculating Ebi, CEbi etc\n")
  #browser()
  nmod <- length(bic)
  model.fits1 <- model.fits
  refit.models <- TRUE
  refit.models <- FALSE
  if (refit.models) {
    model.fits <- as.list(rep(0, nmod))
    for (k in (1:nmod)) {
      if (sum(which[k, ]) != 0) {
        model.fits[[k]] <- ls.print(lsfit(x[, which[k, 
                                                    ], drop = FALSE], y, wt = wt), print.it = FALSE)$coef.table[[1]]
      }
      else model.fits[[k]] <- ls.print(lsfit(rep(1, length(y)), 
                                             y, wt = wt, int = FALSE), print.it = FALSE)$coef.table[[1]]
    }
  }
  Ebi <- rep(0, (nvar + 1))
  SDbi <- rep(0, (nvar + 1))
  CEbi <- Ebi
  CSDbi <- SDbi
  EbiMk <- matrix(rep(0, nmod * (nvar + 1)), nrow = nmod)
  sebiMk <- matrix(rep(0, nmod * (nvar + 1)), nrow = nmod)
  for (i in 1:(nvar + 1)) {
    if ((i == 1) || (sum(which[, (i - 1)] != 0))) {
      for (k in (1:nmod)) {
        if ((i == 1) || (which[k, (i - 1)] == TRUE)) {
          if (i == 1) 
            pos <- 1
          else pos <- 1 + sum(which[k, (1:(i - 1))])
          EbiMk[k, i] <- model.fits[[k]][pos, 1]
          sebiMk[k, i] <- model.fits[[k]][pos, 2]
        }
      }
      Ebi[i] <- as.numeric(sum(postprob * EbiMk[, i]))
      SDbi[i] <- sqrt(postprob %*% (sebiMk[, i]^2) + postprob %*% 
                      ((EbiMk[, i] - Ebi[i])^2))
      if (i == 1) {
        CEbi[i] <- Ebi[i]
        CSDbi[i] <- SDbi[i]
      }
      else {
        sel <- which[, i - 1]
        cpp <- postprob[sel]/sum(postprob[sel])
        CEbi[i] <- as.numeric(sum(cpp * EbiMk[sel, i]))
        CSDbi[i] <- sqrt(cpp %*% (sebiMk[sel, i]^2) + 
                         cpp %*% ((EbiMk[sel, i] - CEbi[i])^2))
      }
    }
  }
  cat("after calculating Ebi, CEbi etc\n")
  #browser()
  dimnames(which) <- list(NULL, dimnames(x)[[2]])
  dimnames(EbiMk) <- dimnames(sebiMk) <- list(NULL, c("Int", 
                                                      dimnames(x)[[2]]))
  if (eval.markers) {
    if (is.null(neval) || neval == 0) 
      neval <- nmod
    d <- matrix(NA, nr = neval, nc = ncol(x))
    for (ii in 1:neval) {
      if (size[ii] > 1) {
        mui <- model.fits[[ii]][1, 1]
        ypi <- as.matrix(x)[, which[ii, ], drop = F] %*% 
          model.fits[[ii]][-1, 1]
        for (jj in 1:ncol(x)) {
          d[ii, jj] <- lm(ypi ~ x[, jj])$coef[2]
        }
      }
      else d[ii, ] <- rep(0, ncol(d))
    }
    d.av <- c(postprob %*% d)
  }
  cat("after eval.markers\n")
  max.size <- max(size)
  mprob.size <- tapply(c(postprob, rep(0, max.size)), c(size, 
                                                        1:max.size) - 1, sum)
  if (is.matrix(probne0)) {
    marker.names <- dimnames(probne0)[[1]]
    probne0 <- c(probne0)
    names(probne0) <- marker.names
  }
  cat("after calculating marginal probabilities for model size\n")
  result <- list(postprob = postprob, r2 = r2, bic = bic, size = (size - 
                                                 1), which = which, probne0 = c(probne0), n = n, postprob.size = mprob.size, 
                 postmean = Ebi, postsd = SDbi, condpostmean = CEbi, condpostsd = CSDbi, 
                 ols = EbiMk, mle = EbiMk, se = sebiMk, n.models = length(postprob), 
                 n.vars = length(probne0))
  if (eval.markers) 
    result$d.av <- d.av
  if (length(prior) == 1 || sqrt(var(prior)) < 1e-06 * mean(prior)) 
    attr(result, "prior") <- prior[1]
  else attr(result, "prior") <- prior
  attr(result, "intercept") <- intercept
  attr(result, "num.imputations") <- num.imputations
  attr(result, "p.sg") <- p.sg
  attr(result, "delta") <- delta
  class(result) <- c("bicreg.qtl", "bicreg")
  result
}



DS.gamma <- function(p){1 + qnorm(1-p/2)*dnorm(qnorm(p/2))/(p/2)}

recalc.bicprobs <- function(obj,old.delta=attr(obj,"delta"),delta=1,n=obj$n,
                            old.prior=attr(obj,"prior"),prior=0.5,p.sg=1){
# calculate estimates of posterior probabilities from BIC
# obj : object returned from biceg.qtl()
# old.delta, delta: old and new values of correction factor delta
# n : length of data (beforre multiple imputations)
# old.prior, prior : old and new values of prior probability for a variable to be in model
#                    (must be scalars)
# p.sg : proportion sampled in the tails
#browser()
   if(length(old.prior) >1 || length(prior) > 1)
     stop("vector prior not supported")
   if(p.sg <0 || p.sg >1)stop("p.sg not a valid proportion")
   kmax <- ncol(obj$which)
   nvars <- length(obj$probne0)
   old.prior <- dbinom(0:kmax,size=kmax,prob=old.prior)
   new.prior <- dbinom(0:kmax,size=kmax,prob=prior)
   gammap <- DS.gamma(p.sg)
   bic.delta <- obj$bic + (delta-old.delta)*obj$size*log(n * gammap)
   bic.delta.adj <- bic.delta -min(bic.delta)
   postprob.delta <- exp(-bic.delta.adj/2)
   postprob.delta <- postprob.delta/sum(postprob.delta)
   postprob.size <-  tapply(c(rep(0,kmax+1),postprob.delta),
        c(0:kmax,obj$size),sum)
   postprob.size <-  postprob.size* new.prior/old.prior
   postprob.size <-  postprob.size/sum(postprob.size)
   postprob <- postprob.delta *new.prior[obj$size+1]/old.prior[obj$size+1]
   postprob <- postprob/sum(postprob)
#   probne0 <- round(100 * t(obj$which) %*% as.matrix(postprob.delta), 1)
   probne0 <- round(100 * t(obj$which) %*% as.matrix(postprob), 1)
   if(!is.null(dimnames(obj$which)[[2]]))names(probne0) <- dimnames(obj$which)[[2]]
#browser()
   list(delta=delta,prior=prior, bic=bic.delta,postprob=postprob.delta,probne0=c(probne0),
        postprob.size=postprob.size)
}

#summary.bicreg.qtl <- function(object,print.it=TRUE){}
summary.bicreg.qtl <- function(object,...){
  nmodels <- length(object$r2)
  # mimum marginal probability to print
  min.marker.prob <- 0.05
  args <- list(...)
  if(!is.null(args$nbest))nmodels <- min(nmodels,args$nbest)
  if(!is.null(args$min.marker.prob))min.marker.prob <- args$min.marker.prob
  model.stats <-cbind(1*(object$which==TRUE), R2=object$r2,BIC=object$bic,
                      postprob=object$postprob,cumprob=cumsum(object$postprob))[1:nmodels,]
  dimnames(model.stats)[[1]] <- 1:nrow(model.stats)

#  cat("R-squared, BIC, and approximate posterior probabilties for individual models:\n")
#  print(model.stats)
#  browser()
  variable.columns <- 1:ncol(object$which)
  variable.mprobs <- sapply(as.list(variable.columns),function(i,m){sum(m[,i]*m[,"postprob"])},
       m=model.stats)
  names(variable.mprobs) <- dimnames(object$which)[[2]]
#  cat("marginal probabilities for individual variables\n")
#  print(variable.mprobs)
  summ <- list(model.stats=model.stats,
               size.marginal.probabilities=object$postprob.size,
               variable.marginal.probabilities=variable.mprobs[variable.mprobs > min.marker.prob])
  attr(summ,"prior") <- attr(object,"prior")
  attr(summ,"intercept") <- attr(object,"intercept")
  attr(summ,"strict") <- attr(object,"strict")
  attr(summ,"OR") <- attr(object,"OR")
  class(summ) <- "summary.bicreg.qtl"
  summ
}


print.summary.bicreg.qtl <- function(x,...){
  cat("R-squared, BIC, and approximate posterior probabilties for individual models:\n")
  print(x$model.stats)
  cat("marginal probabilities for model sizes\n")
  print(x$size.marginal.probabilities)
  cat("marginal probabilities for individual variables\n")
  print(x$variable.marginal.probabilities)
#browser()
  printAttributes(x)
  invisible(x)
}


printAttributes <- function(obj,exclude.attributes=c("names","class")){
  attrs <- attributes(obj)
  for(ii in seq(along=exclude.attributes)){
    attrs[[exclude.attributes[ii]]] <- NULL
  }
  lapply(as.list(seq(along=attrs)), function(ii, attrs){
     at1 <- attrs[[ii]]
     nam1 <- names(attrs)[[ii]]
     cat(paste("attr,\"",nam1,"\")\n",sep=""));print(at1)
  }, attrs=attrs)
  invisible(attrs)
} 




map.function.rec <- function(d,method="haldane"){
 switch(method,
        "haldane"=1/2*(1-exp(-2*d)),
        "kosambi"=(1-exp(-4*d))/(2*(1+exp(-4*d))))
}



map.function.dist <- function(r,method="haldane"){
 switch(method,"haldane"=-1/2*log(1-2*r),
  "kosambi"=1/4*log((1+2*r)/(1-2*r)))
}


scan.string <- function (str, what = numeric(), sep = "", ...) 
{
    if (is.null(str)) 
        return(NULL)
    file <- tempfile("tmpstr")
    on.exit(system(paste("rm", file)))
    cat(paste(str, collapse = sep), file = file)
    scan(file, what = what, sep = sep, ...)
}
