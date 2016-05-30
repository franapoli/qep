
dim.qepmix <- function(qmix)
    return(c(nrow.qepmix(qmix),ncol.qepmix(qmix)))
ncol.qepmix <- function(qmix)
    return(length(dimnames(qmix)[[2]]))
nrow.qepmix <- function(qmix)
    return(length(dimnames(qmix)[[1]]))

qepmix <- function(...)
{
    qdatal <- list(...)
    if(!all(sapply(qdatal, is.qep)))
        stop("all items in qdatal must be of class qep.")

    dict <- unique(do.call(c, lapply(qdatal, rownames)))
    dict <- dict[!is.na(dict)]
    
    return(structure(qdatal, obs=dict, class="qepmix"))
}


`[.qepmix` <- function(qm, v, subslice=T, na.rm=T, collapse=T)
 {
     if(!subslice) {
         v2 <- lapply(v, function(x) rel2abs(qm, x, 1:ncol(qm[[x]])))
         if(length(v)>1)
             v <- do.call(c,v2)
         else v <- v2
     }

     idx <- abs2rel(qm, v)
     qeps <- unique(idx[,1])

     k <- 1

     if(na.rm) {
         nonNAset <- rownames(qm[[qeps[1]]])
         if(length(qeps)>1)
             for(i in 2:length(qeps))
                 nonNAset <- intersect(nonNAset, rownames(qm[[qeps[i]]]))
     } else nonNAset <- attr(qm, "obs")
     
     resmat <- list()
     for(i in 1:length(qeps))
     {
         q <- qeps[i]

         if(na.rm) {
             sub <- qm[[q]][nonNAset, idx[idx[,1]==q,2], drop=F]
             resmat[[i]] <- qep(sub)
         } else
             {
                 sub <- qm[[q]][, idx[idx[,1]==q,2], drop=F]
                 resmat[[i]] <- matrix(NA,length(nonNAset),ncol(sub))
                 rownames(resmat[[i]]) <- nonNAset
                 colnames(resmat[[i]]) <- colnames(sub)
                 resmat[[i]][rownames(sub),] <- sub
             }         
     }
          
     if(collapse)
         if(na.rm)
             return(qep(do.call(cbind, resmat))) else {
                 return(do.call(cbind, resmat))
             }
     
     if(na.rm)
         return(do.call(qepmix, resmat)) else return(resmat)
 }


dims <- function(qm)
{
    dm <- sapply(qm, dim)
    colnames(dm) <- 1:ncol(dm)
    rownames(dm) <- c("Obs.","Cond.")
    return(dm)
}

"dimnames<-.qepmix" <- function(qm, value)
{
    len <- sum(sapply(qm, ncol))
    if(len!=length(value))
        stop("qm number of samples and nms length differ")
    k <- 1
    for(i in 1:length(qm)) {
        n <- ncol(qm[[i]])
        colnames(qm[[i]]) <- value[k:(k+n-1)]
        k <- k+n
    }
    return(qm)
}

"dimnames.qepmix" <- function(qm, i=NULL)
{
    return(list(
        observations=attr(qm, "obs"),
        conditions=do.call(c,sapply(qm, colnames))
        ))
}

abs2rel <- function(qmix, idx)
{
     idx2d <- function(breaks, v)
     {
         idx <- matrix(NA,length(v),2)
         for(i in 1:length(v)) {
         x <- sum(v[i]>breaks)+1
         if(x>1)
             y <- v[i]-breaks[x-1] else y <- v[i]
         idx[i,] <- c(x,y)
         }
         return(idx)
     }

     breaks <- cumsum(sapply(qmix, ncol))
     idx <- idx2d(breaks, idx)
     return(idx)
}


rel2abs <- function(qmix, i, v=NULL)
{
    if(is.null(v))
        v <- 1:ncol(qmix[[i]])
    l <- cumsum(c(0, sapply(qmix, ncol)))
    return(v+l[i])
}


dist.qepmix <- function(qmix, ...)
{
    distmat <- matrix(NA, ncol(qmix), ncol(qmix))
    rownames(distmat) <- colnames(distmat) <- colnames(qmix)
    l <- cumsum(c(0,sapply(qmix, ncol)))
    for(i in 1:length(qmix)) {
        idx <- (l[i]+1):(l[i+1])
        distmat[idx,idx] <- as.matrix(dist.qep(qmix[[i]], ...))
    }

    for(i in 1:(length(qmix)-1))
        for(j in (i+1):length(qmix)) {
            umix <- qmix[c(i, j), subslice=F, collapse=F]
            subdist <- dist.qep(umix[[1]], umix[[2]], ...)
            distmat[rel2abs(qmix,i), rel2abs(qmix,j)] <- subdist
            distmat[rel2abs(qmix,j), rel2abs(qmix,i)] <- t(subdist)
    }

    return(distmat)
}


qmix.subdists <- function(qmix,d)
  {
    n <- length(qmix)
    bounds <- c(1,sapply(qmix,ncol)+1)
    
    res <- list()
    for(i in 1:n){
      res[[i]] <- list()
      for(j in i:n)
        res[[i]][[j]] <- as.vector(d[bounds[i]:(bounds[i+1]-1), bounds[j]:(bounds[j+1]-1)])
      }

    return(res)
  }

qmix.plotsubdists <- function(d)
  {
    n=length(d)
    par(mfrow=c(n,n))
    for(i in 1:n){
      for(j in 1:n)
        if(j<i)
          plot.new() else {
        hist(d[[i]][[j]], 100)
      }
      }    
  }
