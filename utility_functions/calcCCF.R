# Function to calculate cancer cell fraction for a large data frame of mutations.
# Efficient implementation of code from Mcgranahan et al, Science

calcCCF <- function(muts.all = muts.all, pheno=pheno) {
  # Copy of mutation DF, extracting required fields
  m.tmp <- muts.all[,c("alt.reads", "tumour.reads", "patient", "CNt")]
  m.tmp$purity <- pheno$purity[match(as.character(m.tmp$patient), pheno$Sample.Number..WGS.)]
  # Pick out rows with NAs
  m.tmp$ccf.na <- apply(m.tmp[,c('alt.reads', 'tumour.reads', 'CNt')], 1, function(x) {any(is.na(x))})
  
  # tmp <- m.tmp[which(!m.tmp$ccf.na),][1:100,]
  f.function <- function (c,purity,local.copy.number)
  {
    return(min(c((purity*c) / (2*(1-purity) + purity*local.copy.number),1)))
  }  
  
  # Create this distribution for all unique rows (patient doesn't matter here)
  m.tmp$uid <- paste(m.tmp$alt.reads, m.tmp$tumour.reads, m.tmp$CNt, m.tmp$purity, sep='-')
  tmp.uniq <- m.tmp[which(!is.na(m.tmp$ccf.na)),]
  tmp.uniq <- tmp.uniq[which(!duplicated(m.tmp$uid)),]
  s <- seq(0.01,1,length.out=100)
  dists <- lapply(1:dim(tmp.uniq)[1], function(i) {
    x <- dbinom(tmp.uniq$alt.reads[i],tmp.uniq$tumour.reads[i], prob=sapply(s,f.function,tmp.uniq$purity[i],tmp.uniq$CNt[i]))
    names(x) <- s
    return(x)
  })
  sel <- which(sapply(dists, min) == 0)
  for(i in sel) {
    x <- dists[[i]]
    x[length(x)] <- 1
    dists[[i]] <- x
  }
  names(dists) <- tmp.uniq$uid
  
  prob <- 0.95
  dists.m <- lapply(dists, function(x) {
    xnorm   <- x/sum(x)
    xsort   <- sort(xnorm, decreasing = TRUE)
    xcumLik <- cumsum(xsort)
    n = sum(xcumLik < prob) + 1
    LikThresh <- xsort[n]
    cint  <- x[xnorm >= LikThresh]
    cellu <- as.numeric(names(cint))
    m     <- cellu[which.max(cint)]
    m
  })
  dists.clonal <- lapply(dists, function(x) {
    xnorm   <- x/sum(x)
    xsort   <- sort(xnorm, decreasing = TRUE)
    xcumLik <- cumsum(xsort)
    n = sum(xcumLik < prob) + 1
    LikThresh <- xsort[n]
    cint  <- x[xnorm >= LikThresh]
    cellu <- as.numeric(names(cint))
    l.t   <- cellu[1]
    r.t   <- cellu[length(cellu)]
    (l.t <= 1 & r.t >= 1)
  })
  
  # Match back to m.tmp
  m.tmp$ccf <- as.numeric(dists.m[m.tmp$uid])
  m.tmp$is.clonal <- as.logical(dists.clonal[m.tmp$uid])
  
  muts.all$ccf <- m.tmp$ccf
  muts.all$is.clonal <- m.tmp$is.clonal
  
  return(muts.all)
}