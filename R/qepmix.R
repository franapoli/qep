
dim.qepmix <- function(qmix)
    return(c(nrow(qmix),ncol(qmix)))
ncol.qepmix <- function(qmix)
    return(length(dimnames(qmix)[[2]]))
nrow.qepmix <- function(qmix)
    return(length(dimnames(qmix)[[1]]))

qepmix <- function(qdatal)
{
    if(!all(sapply(qdatal, is.qep)))
        stop("all items in qdatal must be of class qep.")

    dict <- sort(unique(do.call(c, lapply(qdatal, rownames))))
    dict <- dict[!is.na(dict)]
    
    return(structure(qdatal, genes=dict, class="qepmix"))
}


unifmix <- function(qm, v, subslice=F, na.rm=T)
 {
     if(!subslice) {
         v2 <- sapply(v, function(x) rel2abs(qmix, x, 1:ncol(qmix[[x]])))
         if(length(v)>1)
             v <- do.call(c,v2)
         else v <- v2
     }

     idx <- abs2rel(qm, v)
     qeps <- unique(idx[,1])

     k <- 1

     if(na.rm) {
         nonNAset <- rownames(qmix[[qeps[1]]])
         if(length(qeps)>1)
             for(i in 2:length(qeps))
                 nonNAset <- intersect(nonNAset, rownames(qm[[qeps[i]]]))
     } else nonNAset <- attr(qm, "genes")

     resmat <- list()
     for(i in 1:length(qeps))
     {
         q <- qeps[i]

         if(na.rm) {
             sub <- qm[[q]][nonNAset, idx[idx[,1]==q,2], drop=F]
             } else sub <- qm[[q]][, idx[idx[,1]==q,2], drop=F]

         resmat[[i]] <- qep(sub)
     }
     
     return(qepmix(resmat))
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
        genes=attr(qm, "genes"),
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
    l <- cumsum(c(0,sapply(qmix, ncol)))
    for(i in 1:length(qmix)) {
        idx <- (l[i]+1):(l[i+1])
        distmat[idx,idx] <- as.matrix(dist(qmix[[i]], method="bsf"))
    }

    for(i in 1:(length(qmix)-1))
        for(j in 2:length(qmix)) {
            umix <- unifmix(qmix, c(i, j))
            subdist <- as.matrix(dist(umix[[1]], umix[[2]], "bsf"))
            distmat[rel2abs(qmix,i), rel2abs(qmix,j)] <- subdist
            distmat[rel2abs(qmix,j), rel2abs(qmix,i)] <- t(subdist)
    }

    return(distmat)
}
