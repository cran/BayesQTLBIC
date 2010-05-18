Unique.Rows <- function(x){
  x <- as.matrix(x)
  (1:nrow(x))[!duplicated(apply(x,1,function(u)paste(u,collapse="\r")))]
}
