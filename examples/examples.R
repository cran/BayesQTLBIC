#l <- dget("/home/rod/clients/genomics/bayesqtl/one-qtl_bc1.dat")
#bc200 <- ex1.df <- data.frame(y=l$y,m=l$m)
#bicreg.qtl(y=l$y,x=l$m,OR=1000,nbest=10,nvmax=4,prior=0.05)

library(leaps)


# back-cross progeny with QTL at positions 40,50,80 cM on chromosome 1
# and 30,55 cM on chromosome 2
ex1.marker.pos <- seq(5,105,by=10)
set.seed(1234)
ex1.qtldata <- sim.bc.progeny(n=1200,Vp=c(0.1,0.2,0.3,0.15,0.25)/2,
      marker.pos=list(chrom=rep(1:2,rep(length(ex1.marker.pos),2)),
      pos=rep(ex1.marker.pos,2)),qtl.pos=list(chrom=rep(1:2,c(3,2)),
                                pos=c(40,50,80,30,55)))
str(ex1.qtldata)
chrom <- rep(1:2,rep(length(ex1.marker.pos),2))
# separate analyses by chromosome
# n=1200, Chromosome 1: 
ex1n1200c1.bicreg <- bicreg.qtl(x=ex1.qtldata$x[,chrom==1],y=ex1.qtldata$y,OR=1000,
                                nbest=10,nvmax=5,prior=0.2,keep.size=1)
# 23 models account for 99 % of the probability
cumsum(ex1n1200c1.bicreg$postprob)
# 2 QTL in coupling at 40,50cM can't be resolved
summary(ex1n1200c1.bicreg,nbest=23)

recalc.bicprobs(ex1n1200c1.bicreg,old.prior=0.2,prior=0.1)


# 75% probability of selecting true model to within resolution of marker map
w1 <- ex1n1200c1.bicreg$which
is.true.modelc1 <- function(x){
  x <- 1*x
  all(x[c(1,2,3,7,10:11)]==0) &&  
  (x[4]==1 || x[5]==1) && (x[5]==1 || x[6]==1) && (x[8]==1 || x[9]==1)  }
apply(w1,1,is.true.modelc1)
sum(apply(w1,1,is.true.model1)*ex1n1200c1.bicreg$postprob)

# n=1200, Chromosome 2: 
ex1n1200c2.bicreg <- bicreg.qtl(x=ex1.qtldata$x[,chrom==2],y=ex1.qtldata$y,OR=1000,
                                nbest=10,nvmax=5,prior=0.2,keep.size=1)
# 11 models account for 99 % of the probability
cumsum(ex1n1200c2.bicreg$postprob)
# QTL at m6 strongly detected
# QTL between m3 and m4 represented by models with one or more
# of m3 and m4 selected
summary(ex1n1200c2.bicreg,nbest=12)
# 94% probability of selecting true model to within resolution of marker map
w2 <- ex1n1200c2.bicreg$which
is.true.modelc2 <- function(x){
  x2 <- 1*x
  all(x2[c(1,2,5,7:11)]==0) &&  x2[6]==1 &&
  (x2[3]==1 || x2[4]==1)   }
sum(apply(w2,1,is.true.modelc2)*ex1n1200c2.bicreg$postprob)

# combined analysis (brute force)
ex1n1200.bicreg <- bicreg.qtl(x=ex1.qtldata$x,y=ex1.qtldata$y,OR=1000,
                              nbest=10,nvmax=5,prior=0.2,keep.size=1)
ex1n1200.bicreg$postprob
cumsum(ex1n1200.bicreg$postprob)
# 7 models explain 99% of the probability
summary(ex1n1200.bicreg,nbest=12,min.marker.prob=0.01)

# same example with ~5% missing values, multiple imputation
ex1.qtldata$x.mv <- ex1.qtldata$x
num.missing <- round(0.05*length(ex1.qtldata$x))
mv.pos <- sample(1:length(c(ex1.qtldata$x)),size=num.missing,replace=F)
ex1.qtldata$x.mv[mv.pos] <- NA

#
if(F){
na.rows <- apply(ex1.qtldata$x.mv,1,function(u)any(is.na(u)))
(1:nrow(ex1.qtldata$x.mv))[na.rows]
ex1.qtldata$x.mv[5:7,]
ex1.qtldata$x.mv[399:413,]
x15 <- impute.marker.matrix(ex1.qtldata$x.mv[1:7,1:9],d=diff(ex1.marker.pos),num.imputations=15)
x300 <- impute.marker.matrix(ex1.qtldata$x.mv[5:7,1:8],d=diff(ex1.marker.pos),num.imputations=300)
apply(matrix(x300[,3],nr=3),1,table)
}

x.imp10 <- impute.marker.matrix(ex1.qtldata$x.mv,d=diff(ex1.marker.pos),num.imputations=10)

ex1n1200.imp10.bicreg <- bicreg.qtl(x=x.imp10,y=rep(ex1.qtldata$y,10),OR=1000,num.imputations=10,
                              nbest=10,nvmax=5,prior=0.2,keep.size=1)
summary(ex1n1200.imp10.bicreg,nbest=12,min.marker.prob=0.01)

w2i10 <-  ex1n1200.imp10.bicreg$which

is.true.model2i10 <- function(x){
  x <- 1*x
  all(x[c(1,2,5,7:11)]==0) &&  x[6]==1 &&
  (x[3]==1 || x[4]==1)   }
sum(apply(w2i10,1,is.true.model2i10)*ex1n1200c2.imp10.bicreg$postprob)


# same example with virtual markers, multiple imputation




# same example with  missing values, n=300

# same example with virtual marker



# a multi-chromosome example
# calculate marginal probabilities
# example with n=300 from QTL co-location paper
load("/home/rod/clients/genomics/QTL+assoc_mapping/example3/ex3n300a.data.rda")
str(ex3n300a.data)
# get x (markers) and y (trait) from QTL Cartographer data.

m <- matrix(scan.string(names(ex3n300a.data$Markers),sep="m",what=""),nr=2)
scan.string(m[1,],sep="c")
chrom <- as.numeric(substring(m[1,],2,nchar(m[1,])))
marker <- as.numeric(m[2,])
x <- sapply(ex3n300a.data$Markers,c)
y <- ex3n300a.data$Trait$t1
nchrom <- length(sort(chrom.levels <- unique(chrom)))

chrom.fits <- list()
for(ii in seq(along=chrom.levels)){
  cat(paste("*** chromosome",ii,"***","\n"))
  ci <- chrom.levels[ii]
  chrom.sel <- chrom==ci
  chrom.fits[[ii]] <- bicreg.qtl(x[,chrom.sel],y, prior=0.1,nbest=20,nvmax=3)
}
class(chrom.fits) <- "bicreg.mqtl"

summary(chrom.fits[[1]], nbest=20)
summary(chrom.fits[[2]], nbest=20)

### sample randomly and independently from models for each other chromosome
#nmodels <- sapply(chrom.fits,function(u)length(u$r2))
set.seed(12345)
nsim <- 20
mWhich <- lapply(chrom.fits, function(fitii){
#browser()
  pprobii <- fitii$postprob
  nmodii <- length(pprobii)
  fitii$which[sample(1:nmodii,size=nsim,prob=pprobii,replace=TRUE),]})
mWhich
nchrom <- length(chrom.fits)
u <- sapply(mWhich, function(m)c(t(m)))
length(u)/(nchrom*nsim)
u1 <- array(c(u),dim=c(length(u)/(nchrom*nsim),nsim,nchrom))
u2 <- t(matrix(aperm(u1,c(1,3,2)),ncol=nsim))
mWhich <- u2


if(F){
dim(u2)
for(ii in 1:nchrom){
 print(any(u2[,16*(ii-1)+1:16] != mWhich[[ii]]))
}
any(u2[,1:16] != mWhich[[1]])
any(u2[,2*16+1:16] != mWhich[[3]])
any(u2[,11*16+1:16] != mWhich[[12]])
}


mWhich <- sample.bicreg.qtl.models(chrom.fits,nsim=50)
summary(mres,nbest=38)
class(mres)
plot(mres)
mres50 <- mres
class(mres) <- "bicreg"
plot(mres)

mWhich200 <- sample.bicreg.qtl.models(chrom.fits,nsim=200)
mres200 <- bicreg.models(x=x,y=y,which=mWhich200,prior=0.1)

summary(mres200,nbest=38,min.marker.prob=0.05)
# marginal probabilities for markers greater than 10%
mres200$probne0[mres200$probne0 > 5]


# 33 models explain 99% of the probability from this sample.
mres200$cumprob[1:38]



test.mres <- bicreg.models(x=x,y=y,which=mWhich,prior=0.1)


mWhich1000 <- sample.bicreg.qtl.models(chrom.fits,nsim=1000)
mres1000 <- bicreg.models(x=x,y=y,which=mWhich1000,prior=0.1)

# 33 models explain 99% of the probability, still
summary(mres1000,nbest=10,min.marker.prob=0.05)
# timeseries (genome) plot of marginal probabilities for markers
xp <- chrom + (marker -1 )/16
plot(xp,mres1000$probne0,type="l",xaxt="n",
     xlab="genome position", ylab="marginal probability per marker (*100)")
axis(1,at=0.5+1:12,lab=paste("c",1:12),tck=par()$tck*1.2)
axis(1,at=xp,lab=rep("",length(xp)))
abline(v=1:13,lty=2)
