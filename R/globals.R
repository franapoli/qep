    
ncol <- function(object, ...) UseMethod("ncol")
ncol.default <- base::ncol
nrow <- function(object, ...) UseMethod("nrow")
nrow.default <- base::nrow
dist <- function(object, ...) UseMethod("dist")
if(require(stats))
    dist.default <- stats::dist

schementr <- function(scheme)
{
    if(abs(sum(scheme)-1)>10^-15)
        stop("scheme values must sum to 1")
    -sum(scheme*log(scheme))
}

bsf.dist.row <- function(v1, V2, nbins, SF, maxD)
{
    return(
        apply(V2, 2, function(x) bsf.dist(v1, x, nbins, SF, maxD))
    )
}

jacc.dist.row <- function(v1, V2)
{
    return(
        apply(V2, 2, function(x) sum(abs(v1 - x)))
    )
}

mantra.dist <- function(v1, v2, ng)
    {
        gseadist <- function(v1,v2, binidx)
            {
                n <- length(v1)
                v1inv2 <- v2[match(binidx, v1)]
                binvect <- 1:n %in% v1inv2
                NR <- length(binidx)
                hits <- cumsum(binvect)/NR
                miss <- cumsum(!binvect)/(n-NR)
                w <- which.max(abs(hits-miss))
                return(hits[w]-miss[w])
            }

        n <- length(v1)
        d1 <- gseadist(v1,v2,1:ng)
        d2 <- gseadist(v1,v2,(n-ng+1):n)
        d3 <- gseadist(v2,v1,1:ng)
        d4 <- gseadist(v2,v1,(n-ng+1):n)
        Davg <- ((d1-d2)/2 + (d3-d4)/2)/2
        Dmax <- max((d1-d2)/2, (d3-d4)/2)        
        
        return(1 - Davg)
    }
    
    
mantra.dist.row <- function(v1, V2, ng)
{
    return(
        apply(V2, 2, function(x) mantra.dist(v1, x, ng))
    )    
}

bsf.dist <- function(v1, v2, nbins, SF=NULL, maxD=NULL) {
    if(is.null(maxD))
        maxD <- sum(abs(1:nbins-(nbins:1)))
    if(is.null(SF))
        SF <- abs(row(matrix(NA,nbins,nbins))-col(matrix(NA,nbins,nbins)))    
    
    v1 <- as.factor(v1)
    v2 <- as.factor(v2)
    ints <- table(v1,v2)
    ridx <- as.numeric(rownames(ints))
    cidx <- as.numeric(colnames(ints))    
                           
    unions <- outer(table(v1),table(v2),'+')-ints
    JSF <- SF[ridx,cidx]*(ints/unions)
    return(sum(JSF)/maxD)
}


list2distmat <- function(d, dnames)
  {
    n <- length(d)+1
    m <- matrix(0, n, n)
    for(i in 1:(n-1))
        m[(i+1):n,i] <- d[[i]]
    rownames(m) <- colnames(m) <- dnames
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
