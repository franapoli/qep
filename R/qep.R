
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

    if(distfname == "bsf.dist.row") {
        nbins <- max(qep1)
        SF <- abs(row(matrix(NA,nbins,nbins))-col(matrix(NA,nbins,nbins)))
        maxD <- sum(abs(1:nbins-nbins:1))
        distpar <- list(nbins = nbins, SF = SF, maxD = maxD)
    }    
  
    if(is.null(qep2)) {
        dists <- matrix(NA, ncol(qep1), ncol(qep1))
        rownames(dists) <- colnames(dists) <- colnames(qep1)
        nsteps <- ncol(qep1)-1
        
        expr <- expression({
            if(!is.null(callbackf))
                do.call(callbackf, list(i))
            idx <- (i+1):ncol(qep1)
            do.call(distf, c(list(v1 = qep1[,i], V2 = qep1[,idx,drop=F]),
                             distpar))
        })
        
        if(parallel) {
            dists <- foreach(i=1:nsteps, .export=distfname) %dopar% eval(expr)
        } else dists <- foreach(i=1:nsteps) %do% eval(expr)

        return(list2distmat(dists))

  
    } else {
        warning("implement me")
                        ## dists <- matrix(NA, ncol(qep1), ncol(qep2))
                        ## rownames(dists) <- colnames(qep1)
                        ## colnames(dists) <- colnames(qep2)
                        ## for(i in 1:ncol(qep1))
                        ##     dists[i,] <- apply(qep2,2,bsf,v1=qep1[,i])
                        ## return(dists)
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


