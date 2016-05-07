\name{logistic.ep}
\alias{logistic.ep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
logistic.ep(nbins, steep, width, baseline, zero, trim = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nbins}{
%%     ~~Describe \code{nbins} here~~
}
  \item{steep}{
%%     ~~Describe \code{steep} here~~
}
  \item{width}{
%%     ~~Describe \code{width} here~~
}
  \item{baseline}{
%%     ~~Describe \code{baseline} here~~
}
  \item{zero}{
%%     ~~Describe \code{zero} here~~
}
  \item{trim}{
%%     ~~Describe \code{trim} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (nbins, steep, width, baseline, zero, trim = T) 
{
    l1 <- steep
    hsh <- width
    vsh <- baseline
    rot <- zero
    odd <- nbins\%\%2
    nbins <- ceiling(nbins/2)
    f <- plogis(seq(-l1, l1, length = nbins) + hsh) + vsh
    f <- (1 - vsh)/(max(f) - min(f)) * (f - max(f)) + 1
    f2 <- f[length(f):1]
    if (odd) 
        f <- c(f[1:(length(f) - 1)], f2)
    else f <- c(f, f2)
    shiftf <- function(f, rot) {
        l <- length(f)
        if (!trim) {
            f <- f[rot:(l + rot - 1)\%\%l + 1]
        }
        else {
            if (rot > 0) {
                f <- c(rep(f[1], rot), f)[1:length(f)]
            }
            else {
                f <- tail(c(f, rep(tail(f, 1), abs(rot))), length(f))
            }
            f <- f/sum(f)
        }
        return(f)
    }
    if (is.null(dim(zero))) {
        f <- shiftf(f, rot)
    }
    else {
        rots <- 0.5 - apply(zero, 2, function(x) sum(x > 0)/length(x))
        f <- sapply(round(rots * nbins), shiftf, f = f)
    }
    return(f)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
