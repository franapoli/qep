\name{quantize.ep}
\alias{quantize.ep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
quantize.ep(data, scheme, decreasing = T, verbose = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{scheme}{
%%     ~~Describe \code{scheme} here~~
}
  \item{decreasing}{
%%     ~~Describe \code{decreasing} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
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
function (data, scheme, decreasing = T, verbose = T) 
{
    qep_checknames(data)
    if (decreasing) {
        data <- apply(-data, 2, rank, ties.method = "random")
    }
    else data <- apply(data, 2, rank, ties.method = "random")
    data <- (data - 1)/(nrow(data) - 1)
    if (!is.null(dim(scheme))) {
        if (abs(sum(scheme) - ncol(scheme)) > (10^-15) * ncol(scheme)) 
            stop("All schemes must sum to 1")
        binned <- matrix(NA, nrow(data), ncol(data))
        for (i in 1:ncol(data)) {
            cs <- cumsum(scheme[, i])
            cs[length(cs)] <- 1
            binned[, i] <- cut(data[, i], c(0, cs), F, T)
        }
    }
    else {
        cs <- cumsum(scheme)
        cs[length(cs)] <- 1
        binned <- apply(data, 2, cut, c(0, cs), F, T)
    }
    rownames(binned) <- rownames(data)
    colnames(binned) <- colnames(data)
    return(qep(binned))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
