\name{dist.qep}
\alias{dist.qep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dist.qep(qep1, qep2 = NULL, method = "manhattan", l = NULL, parallel = F, verbose = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{qep1}{
%%     ~~Describe \code{qep1} here~~
}
  \item{qep2}{
%%     ~~Describe \code{qep2} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{l}{
%%     ~~Describe \code{l} here~~
}
  \item{parallel}{
%%     ~~Describe \code{parallel} here~~
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
function (qep1, qep2 = NULL, method = "manhattan", l = NULL, 
    parallel = F, verbose = T) 
{
    nbins <- max(qep1)
    SF <- abs(row(matrix(NA, nbins, nbins)) - col(matrix(NA, 
        nbins, nbins)))
    maxD <- sum(abs(1:nbins - nbins:1))
    maxD <- 1
    if (is.null(qep2)) {
        d <- switch(method, manhattan = stats::dist(qep1, "manhattan"), 
            bsf = {
                dists <- matrix(NA, ncol(qep1), ncol(qep1))
                rownames(dists) <- colnames(dists) <- colnames(qep1)
                nsteps <- ncol(qep1) - 1
                dists <- foreach(i = 1:nsteps, .export = "bsf.dist") \%dopar\% 
                  {
                    cat("\r row", i, " of ", nsteps)
                    sapply((i + 1):ncol(qep1), function(x) bsf.dist(qep1[, 
                      i], qep1[, x], nbins, SF, maxD))
                  }
                cat("\n")
                return(list2distmat(dists))
            })
    }
    else {
        d <- switch(method, manhattan = stats::dist(qep1, "manhattan"), 
            bsf = {
                dists <- matrix(NA, ncol(qep1), ncol(qep2))
                rownames(dists) <- colnames(qep1)
                colnames(dists) <- colnames(qep2)
                for (i in 1:ncol(qep1)) dists[i, ] <- apply(qep2, 
                  2, bsf, v1 = qep1[, i])
                return(dists)
            })
    }
    return(d)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
