
shiftf <- function(f, zero, trim=T) {

    doit <- function(f, rot, trim) {
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
        f <- doit(f,zero,trim)
    } else {
      f <- sapply(round(zero*length(f)), shiftf, f=f, trim=trim)
    }
    
    return(f)
}


unif.qs <- function(nbins)
{
    uni <- rep(1/nbins, nbins)
    return(uni)
}


norm.qs <- function(nbins, sd, baseline, zero=0, trim=T)
{
    vsh <- baseline
    uni <- seq(0,1,length=nbins)
    f <- dnorm(seq(-1,1,length=nbins), sd=sd)
    if(max(f)-min(f)>0) {
      f <- (1-vsh)/(max(f)-min(f))*(f-max(f))+1 # range normalization
    } else f <- f/sum(f)
    f <- shiftf(f,zero,trim)
    return(f)
}

logis.qs <- function(nbins, steep, width, baseline, zero=0, trim=T)
{
  if(nbins==1) return(1)
  
    l1 <-  30*steep;
    hsh <- 30*width;
    vsh <- baseline; rot <- zero
  
    odd <- nbins%%2
    nbins <- ceiling(nbins/2)
    f <- plogis(seq(-l1,l1,length=nbins)+hsh)+vsh;

    if(max(f)-min(f)>0) {
      f <- (1-vsh)/(max(f)-min(f))*(f-max(f))+1 # range normalization
    } else f <- f/sum(f)

    f2 <- f[length(f):1]

  if(odd) {
      f <- c(f[1:(length(f)-1)],f2)
  } else f <- c(f,f2)

  f <- shiftf(f,zero,trim)

    return(f)
}
                          
quantize <- function(data, qscheme,
                     decreasing=T, verbose=T)
{
    
    if(decreasing) {
        data <- apply(-data,2,rank,ties.method="random")
        } else data <- apply(data,2,rank,ties.method="random")

    data <- (data-1)/(nrow(data)-1)


    if(!is.null(dim(qscheme))) {
        if(abs(sum(qscheme)-ncol(qscheme))>(10^-15)*ncol(qscheme))
            stop("All schemes must sum to 1")
        if(ncol(qscheme)!=ncol(data))
            stop("data and scheme must have same number of columns")
        binned <- matrix(NA,nrow(data),ncol(data))
        for(i in 1:ncol(data)) {
          cs <- cumsum(qscheme[,i])
          ## make sure 1 is 1 and not 1-10^16...    
          cs[length(cs)] <- 1
          binned[,i] <- cut(data[,i], c(0,cs), F, T)
        }
    } else {
      cs <- cumsum(qscheme)
      ## make sure 1 is 1 and not 1-10^16...    
      cs[length(cs)] <- 1
      binned <- apply(data, 2, cut, c(0,cs), F, T)
    }

    rownames(binned) <- rownames(data)
    colnames(binned) <- colnames(data)
    
    return(qep(binned))
}

