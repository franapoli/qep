
logistic.ep <- function(nbins, steep, width, baseline, zero, trim=T)
{
    l1 <-  30*steep;
    hsh <- 30*width;
    vsh <- baseline; rot <- zero
  
    odd <- nbins%%2
    nbins <- ceiling(nbins/2)
    f <- plogis(seq(-l1,l1,length=nbins)+hsh)+vsh;
    f <- (1-vsh)/(max(f)-min(f))*(f-max(f))+1 # range normalization

    f2 <- f[length(f):1]

    if(odd)
        f <- c(f[1:(length(f)-1)],f2) else f <- c(f,f2)

    shiftf <- function(f, rot) {
        
        l <- length(f)
        if(!trim) {
            f <- f[rot:(l+rot-1) %% l + 1]
        } else {
            if(rot>0) {
                f <- c(rep(f[1], rot), f)[1:length(f)]
            } else {
                f <- tail(c(f, rep(tail(f,1), abs(rot))), length(f))
            }
            f <- f/sum(f)
        }
        return(f)
    }

    if(length(zero)==1)
        {
            f <- shiftf(f,rot)
        } else {
            f <- sapply(round(rot*nbins), shiftf, f=f)
        }

    return(f)
}
                          
quantize.ep <- function(data, scheme,
                        decreasing=T, verbose=T)
{
    qep_checknames(data)
    
    if(decreasing) {
        data <- apply(-data,2,rank,ties.method="random")
        } else data <- apply(data,2,rank,ties.method="random")

    data <- (data-1)/(nrow(data)-1)


    if(!is.null(dim(scheme))) {
        if(abs(sum(scheme)-ncol(scheme))>(10^-15)*ncol(scheme))
            stop("All schemes must sum to 1")
        if(ncol(scheme)!=ncol(data))
            stop("data and scheme must have same number of columns")
        binned <- matrix(NA,nrow(data),ncol(data))
        for(i in 1:ncol(data)) {
          cs <- cumsum(scheme[,i])
          ## make sure 1 is 1 and not 1-10^16...    
          cs[length(cs)] <- 1
          binned[,i] <- cut(data[,i], c(0,cs), F, T)
        }
    } else {
      cs <- cumsum(scheme)
      ## make sure 1 is 1 and not 1-10^16...    
      cs[length(cs)] <- 1
      binned <- apply(data, 2, cut, c(0,cs), F, T)
    }

    rownames(binned) <- rownames(data)
    colnames(binned) <- colnames(data)
    
    return(qep(binned))
}

