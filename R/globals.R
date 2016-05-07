
ncol <- function(object, ...) UseMethod("ncol")
ncol.default <- base::ncol
nrow <- function(object, ...) UseMethod("nrow")
nrow.default <- base::nrow
dist <- function(object, ...) UseMethod("dist")
if(require(stats))
    dist.default <- stats::dist


bsf.dist <- function(v1, v2, nbins=max(c(v1,v2)), SF, maxD=NULL) {
    if(is.null(maxD))
        maxD <- sum(abs(1:nbins-(nbins:1)))
    v1 <- as.factor(v1)
    v2 <- as.factor(v2)
    ints <- table(v1,v2)
    ## ridx <- as.numeric(rownames(ints))
    ## cidx <- as.numeric(colnames(ints))    
    ridx <- as.numeric(rownames(ints))
    cidx <- as.numeric(colnames(ints))    
                           
    unions <- outer(table(v1),table(v2),'+')-ints
    unions[unions==0] <- 1      
    JSF <- SF[ridx,cidx]*(ints/unions)
    return(sum(JSF, na.rm=T)/maxD)
}


list2distmat <- function(d)
  {
    n <- length(d)+1
    m <- matrix(0, n, n)
    for(i in 1:(n-1))
      m[(i+1):n,i] <- d[[i]]
    return(as.dist(m))
  }

synData <- function(nobsuniv, nobs, nconds, nbins, meansd=.05)
  {
    genome <- paste0("O",1:nobsuniv)

    rawdata <- list()
    qdata <- list()
    for(i in 1:length(nobs))
      {
        nr <- nobs[[i]]; nc <- nconds[[i]]
        m <- matrix(rnorm(nr*nc, 0.5),nr,nc)
        m <- apply(m, 2, function(x) x + rnorm(1, sd=meansd))
        rawdata[[i]] <- m
        
        rownames(rawdata[[i]]) <- sample(genome, nr)
        colnames(rawdata[[i]]) <- paste0("P",i,"_C", 1:nc)
        
        s <- logistic.ep(nbins, steep=.5, width=0, baseline=.1, zero=rawdata[[i]])
        
        qdata[[i]] <- quantize.ep(rawdata[[i]], s)
      }
    
    qmix <- qepmix(qdata)

    return(list(rawdata = rawdata, qep = qdata))
  }