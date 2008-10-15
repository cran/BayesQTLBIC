sim.bc.progeny <-
  function(n,Vp=NULL,map.pos,qtl.pos){
    # simulated one backcross progeny with give qtl and markers
    # for multiple linkage groups map.pos,qtl.pos should be lists
    # map.pos$chrom, map.pos$pos etc.    
    if(!is.list(map.pos)){
      map.pos <- list(chrom=rep(1,length(map.pos)), pos=map.pos)
    }
    if(is.list(qtl.pos)){
      if(is.null(Vp))Vp <- qtl.pos$h2qs*qtl.pos$qtl.effect.signs
    }
    if(!is.list(qtl.pos)){
      qtl.pos <- list(chrom=rep(1,length(qtl.pos)), pos=qtl.pos)
    }
    chrom.u <- sort(unique(c(qtl.pos$chrom,map.pos$chrom)))
    if(sum(Vp)>1)stop("QTL variance > total")
    if(length(Vp)!=length(qtl.pos$pos))stop("incompatible QTL variances and positions")
    res <- list()
    for(ic in seq(along=chrom.u)){
      cc <- chrom.u[ic]
      cat("chromosome",cc,"\n")
      map.pos1 <- map.pos$pos[map.pos$chrom==cc]
      qtl.pos1 <- qtl.pos$pos[qtl.pos$chrom==cc]
      Vp1 <- Vp[qtl.pos$chrom==cc]
      res[[cc]] <- list()
      for(ii in 1:n){
        res[[cc]][[ii]] <- sim.1bc.progeny(Vp1,map.pos1,qtl.pos1)
      }
    }
    ## list of matrices, rows <-> markers on each chromosome
    x.list <- lapply(res,function(c1)sapply(c1,function(u)u$marker.values))
    markers.per.chrom <- unlist(lapply(x.list,function(u)nrow(u)))
    x <- x.list[[1]]
    if(length(x.list)>1){
      for(ii in 2:length(x.list)){
        x <- rbind(x,x.list[[ii]])
      }
    }
    ##    x <- t(sapply(res,function(u)u$marker.values))
    x <- t(x)
    col.chrom <- rep(chrom.u,markers.per.chrom)
    col.chrom
    col.marker <- c(unlist(lapply(split(col.chrom,col.chrom),function(u)seq(along=u))))
    dimnames(x) <- list(NULL,
                        paste("c",col.chrom,"m",col.marker,sep=""))

    nq.list <- sapply(lapply(res,function(c1){
      sapply(c1,function(u)length(u$qtl.effect) )}),
                      mean)

    yq.list <- lapply(res[nq.list >0],function(c1){
      y <- sapply(c1,function(u)u$qtl.effect) 
                                        #      if(is.matrix(y))y <- apply(y,2,sum)
      y
    })
    y <- apply(sapply(yq.list,function(u){
      if(!is.matrix(u))u else apply(u,2,sum)}),1,sum)
    err <- sqrt((1 - sum(Vp)))*rnorm(n)
                                        #browser()
    list(x=x,y=y+err)
  }

sim.1bc.progeny <-  function(Vp,map.pos,qtl.pos){
### simulated one backcross progeny with given qtl and markers
                                        # for one linkage group
    positions <- sort(c(map.pos,qtl.pos))
    signs <- sign(Vp)
    Vp <- abs(Vp)
    a <- 2*sqrt(Vp) * signs
    qtl.effects <- cbind(-a/2,a/2)

    ds <- diff(positions)/100
    rs <- map.function.rec(ds,method="kosambi")
    recs <- numeric(length(rs))
    for(ii in seq(along=rs)){
      recs[ii] <- rbinom(n=1,size=1,prob=rs[ii])
    }
    recs <- c(rbinom(n=1,size=1,prob=0.5),recs)
    recs
    phases <- cumprod(ifelse(recs==0,1,-1))
    phases
    marker.values <- ifelse(phases==1,1,2)
    qtl.ind <- match(qtl.pos,positions)
    qtl.value <- marker.values[qtl.ind]
    qtl.effect <- numeric(0)
    if(length(qtl.ind)>0){
      marker.values <- marker.values[-qtl.ind]
      qtl.effect <- qtl.effects[cbind(1:nrow(qtl.effects),qtl.value)]
    }
                                        #if(length(qtl.pos)==0)browser()
    debug <- FALSE
    if(debug){
      print(list(marker.values=marker.values,qtl.value=qtl.value, 
                 qtl.effect =qtl.effect))

      browser()
    }
    list(marker.values=marker.values,qtl.value=qtl.value, 
         qtl.effect =qtl.effect)
  }


if(F){

sim.bc.progeny <- function(n,Vp,map.pos,qtl.pos){
  ## simulated one backcross progeny with give qtl and markers
  ## for one linkage group
  res <- list()
  for(ii in 1:n){
    res[[ii]] <- sim.1bc.progeny(Vp,map.pos,qtl.pos)


  }
  x <- t(sapply(res,function(u)u$marker.values))
  y <- sapply(res,function(u)u$qtl.effect) 
  err <- (1 - Vp)*rnorm(n)
  list(x=x,y=y+err)
}

if(F){
sim.1bc.progeny <- function(Vp,map.pos,qtl.pos){
# simulate one backcross progeny with given qtl and markers
# for one linkage group
positions <- sort(c(map.pos,qtl.pos))
a <- 2*sqrt(Vp)
qtl.effects <- c(-a/2,a/2)
ds <- diff(positions)/100
rs <- map.function.rec(ds,method="kosambi")
recs <- numeric(length(rs))
for(ii in seq(along=rs)){
  recs[ii] <- rbinom(n=1,size=1,prob=rs[ii])
}
recs <- c(rbinom(n=1,size=1,prob=0.5),recs)
recs
phases <- cumprod(ifelse(recs==0,1,-1))
phases
marker.values <- ifelse(phases==1,1,2)
qtl.ind <- match(qtl.pos,positions)
qtl.value <- marker.values[qtl.ind]

list(marker.values=marker.values[-qtl.ind],qtl.value=qtl.value, qtl.effect =qtl.effects[qtl.value])
}
}
}
