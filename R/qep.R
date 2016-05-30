
#'@import foreach

qep_checknames <- function(qdata)
{    
    gnames <- rownames(qdata)
    if(class(qdata)!="matrix")
        stop("qdata must be an object of class matrix.")
    if(is.null(gnames))
        stop(paste0("Please provide observation names as rownames ",
                    "of qdata or through the gnames parameter."))
    if(any(duplicated(gnames)))
        stop("Observation names must be unique.")
    if(any(is.na(gnames)))
        stop("Observation names can not be NA")
}

qep <- function(qdata, nbins = max(qdata), suppressMinWarn = F)
{
    if(!is.integer(qdata)) {
        if(all(round(qdata)==qdata)) {
            qdata <- matrix(as.integer(qdata), nrow(qdata),
                            dimnames=list(rownames(qdata),colnames(qdata)))
        } else stop("qdata must contain integer values.")
    }

    if(any(is.na(qdata)))
        stop("qdata can not contain NAs")

    if(min(qdata)>1 && !suppressMinWarn)
        warning("Minimum value larger than 1")

    qep_checknames(qdata)
    
    dimnames(qdata) <- list(observations=rownames(qdata), conditions=colnames(qdata))
    
    return(
        structure(
            qdata,
            nbins = nbins,
            class = "qep"
        )
    )
}

is.qep <- function(q)
{
    return(class(q)=="qep")
}


dist.qep <- function(qep1, qep2=NULL, distf=bsf.dist.row,
                     verbose=T, parallel=F, callbackf=NULL, ...)
{
    distpar <- list(...)
    distfname <- deparse(substitute(distf))
  
    if(is.null(qep2)) {
        expr <- expression({
            if(!is.null(callbackf))
                do.call(callbackf, list(i))
            idx <- (i+1):ncol(qep1)
            do.call(distf, c(list(v1 = qep1[,i], V2 = qep1[,idx,drop=F]),
                             distpar))
        })
        
        if(parallel) {
            dists <- foreach(i=1:(ncol(qep1)-1),
                             .export=c(distfname, "callbackf", "qep1", "distpar", "distf")
                             ) %dopar% eval(expr)
        } else dists <- foreach(i=1:(ncol(qep1)-1)) %do% eval(expr)
        
        return(list2distmat(dists, colnames(qep1)))
  
    } else {
        expr <- expression({
            if(!is.null(callbackf))
                do.call(callbackf, list(i))
            do.call(distf, c(list(v1 = qep1[,i], V2 = qep2),
                             distpar))
        })        
        
        if(parallel) {
            dists <- foreach(i=1:ncol(qep1),
                             .export=c(distfname, "callbackf", "qep1", "qep2", "distpar", "distf")
                             )%dopar% eval(expr)
        } else dists <- foreach(i=1:ncol(qep1)) %do% eval(expr)

        m <- do.call(rbind, dists)
        rownames(m) <- colnames(qep1)
        colnames(m) <- colnames(qep2)
        return(m)
    }
}


print.qep <- function(x, ...)
{
    nbins <- attr(x, "nbins")
    attr(x, "class") <- NULL
    attr(x, "nbins") <- NULL
    names(dimnames(x)) <- c("Obs.", "Conditions")
    print(x)
    cat(paste("Num. bins:", nbins, "\n"))
}

## `[.qep` <- function(x, ..., compact=F) {
##     if(compact){
##         p <- list(...)
##         attr(x, "class") <- NULL
##         x <- list(...))
##         return(qep(apply(x, 2, rank, ties.method="random")))
##     }
##     NextMethod()
## }

summary.qep <- function(x, nbins=max(x, na.rm=T))
{
    cat("Number of conditions:", ncol(x), "\n")
    cat("Number of observations:", nrow(x), "\n")
    missg <- unique(apply(x, 2, function(a) sum(is.na(a))))
    cat("\n")
    cat("Unique types of binning:\n")

    tabs <- t(apply(x,2, function(a)
        table(factor(a, levels=1:nbins))
        ))

    utabs <- unique(tabs)

    times <- vector("numeric", nrow(utabs))
    for(i in 1:nrow(utabs))
        times[i] <- sum(apply(tabs,1,function(x) all(x==utabs[i,])))

    res <- t(utabs)
    rownames(res) <- paste(1:nrow(res))
    colnames(res) <- times
    res <- res[,order(as.numeric(colnames(res)), decreasing=T), drop=F] 
    
    print(res)

    return(invisible(res))
}


