cor.z <- function(data.matrix) {
  suppressWarnings(cor.matrix <- cor(data.matrix))
  cor.matrix[is.na(cor.matrix)] <- 0
  cor.matrix[cor.matrix==1] <- 0
  cor.matrix <- z.transform(cor.matrix)
  return(cor.matrix)
}

cor.p <- function(pearson.r,n) {
  t.stat <- pearson.r*sqrt((n-2)/(1-pearson.r^2))
  p = pt(t.stat,n-2)
  return(1-p)
}

cor.threshold <- function(cor.p,n) {
  qval <- 1 - cor.p
  t.stat <- qt(qval,n-2)
  temp <- t.stat^2/(n-2)
  r <- sqrt(temp/(temp+1))
  return(r)
}

z.transform <- function(r) {
  z = 0.5*log((1+r)/(1-r))
  return(z)
}

z.inverse <- function(z) {
  r <- (exp(2*z)-1)/(exp(2*z)+1)
  return(r)
}

cor.jack <- function(data.matrix,njack=1,nresample=NULL) {
  nfeat <- dim(data.matrix)[2]
  nsamples <- dim(data.matrix)[1]
  if (njack==1) {
    cor.matrix.jack <- matrix(0,nfeat,nfeat)
    for (i in 1:nsamples) {
      cor.matrix.jack <- cor.matrix.jack + cor.z(data.matrix[-i,])
      cat(".")
    }
    cor.matrix.jack <- cor.matrix.jack/nsamples
    cor.matrix.tilda <- nsamples*cor.z(data.matrix)-(nsamples-1)*cor.matrix.jack
    cat(".\n")
    cor.matrix.tilda <- z.inverse(cor.matrix.tilda)
    diag(cor.matrix.tilda) <- 1
    return(cor.matrix.tilda)
  }
  else {
    # We only do leave-one-out jackknife first
    # To be added...
    return(0)
  }
}
