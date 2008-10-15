impute.marker.matrix <- function(m,d,num.imputations=1, map.function="haldane",marker.values=1:2){
# m a matrix of marker values
res <- list()
if(!missing(num.imputations)) num.imputations <- as.integer(num.imputations)
for(ii in 1:num.imputations){
  res[[ii]] <- apply(m,1,function(u,d1,map.funct,marker.v){
#    browser()
    impute.marker(u,d=d1,map.function=map.funct,marker.values=marker.v)
  },
                   d1=d,map.funct=map.function,marker.v=marker.values)
}
#browser()
invisible( t(do.call(cbind, res)) )
}

impute.marker <-
function(m, d,  map.function = "haldane", marker.values = 1:2){
# randomly impute missing values where marker is missing
# according to probabilities
# m : vector of marker values (integer values  in marker.values)
# d[i] = distance from m[i] to m[i+1], in cM
# map.function : mapping function to use ("haldane" or
# pos[i] = position of m[i] on its linkage group in cM
# m.int <- c(factor(m, levels=marker.values))
# m.int.values <- as.numeric(levels(factor(m.int)))
# marker.values[m.int]
#browser()
  m.int <- match(m,marker.values)
  m.int.values <- seq(along=marker.values)
  if(any(!is.na(m) & is.na(match(m, marker.values)))) 
    stop("illegal value for marker")
  mf <- switch(map.function,
               haldane = function(d)
               map.function.rec(d, method = "haldane"),
               kosambi = function(d)
               map.function.rec(d, method = "kosambi"))
  if(any(is.na(m))) {
    na.ind <- seq(along = m)[is.na(m)]
    if(length(na.ind == 1)) ii <- na.ind else 
    ii <- sample(na.ind, size = length(na.ind), replace = F)
    ok.ind <- seq(along = m)[!is.na(m)]
    for(jj in ii) {
      # each missing marker
      # get non-missing flanking markers
      # case of no flanking markers
      if(!any(ok.ind)) m[jj] <- sample(marker.values,size=1,prob=c(0.5,0.5)) else {
        # exists flanking markers
        if(any(ok.ind < jj)) l1 <- max(ok.ind[ok.ind < jj]) else l1 <- NULL
        if(any(ok.ind > jj)) r1 <- min(ok.ind[ok.ind > jj]) else r1 <- NULL
        if(!is.null(l1))
          dl <- sum(d[l1:(jj - 1)])
        if(!is.null(r1))
          dr <- sum(d[jj:(r1 - 1)])
        ## flanking marker on the left only
        if(!is.null(l1) & is.null(r1)) {
          rl <- mf(dl/100)
          if(m[l1] == marker.values[2])
            rl <- 1 - rl
          m[jj] <- sample(marker.values, size = 1, prob
                          = c(1 - rl, rl))
          m.int[jj] <- match(m[jj],marker.values)
          ok.ind <- c(jj, ok.ind)
        }
        ## flanking marker on the right only
        if(is.null(l1) & !is.null(r1)) {
          rr <- mf(dr/100)
          if(m[r1] == marker.values[2])
            rr <- 1 - rr
          m[jj] <- sample(marker.values, size = 1, prob
                          = c(1 - rr, rr))
          m.int[jj] <- match(m[jj],marker.values)
          ok.ind <- c(jj, ok.ind)
        }
        ## 2 flanking markers
        if(!is.null(l1) & !is.null(r1)) {
          rl <- mf(dl/100)
          rr <- mf(dr/100)
          tbl <- haldane.tbl(rl, rr)
#browser()
          m.int[jj] <- sample(m.int.values, size = 1, prob = 
                              c(tbl[1, m.int[l1], m.int[r1]], tbl[2, m.int[l1], m.int[r1]])
)
          m[jj] <- marker.values[m.int[jj]]
        }
      }
    } # for(jj in ii)
  } # if (any(is.na(m))
  m
}

haldane.tbl <- function(rl,rr){
# tbl[i,j,k] is probability
# marker is i given flanking marker values j,k
# for i,j,k = 1,2
# flanking markers at distance rl,rr  
   tbl <- array(0,dim=c(2,2,2))
   tbl[1,1,1] <- (1-rl)*(1-rr)
   tbl[2,1,1] <- rl*rr
   tbl[1,1,2] <- (1-rl)*rr
   tbl[2,1,2] <- rl*(1-rr)
   tbl[1,2,1] <- rl*(1-rr)
   tbl[2,2,1] <- (1-rl)*rr
   tbl[1,2,2] <- rl*rr
   tbl[2,2,2] <- (1-rl)*(1-rr)
   tbl
}
