\name{bsf.dist}
\alias{bsf.dist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bsf.dist(v1, v2, nbins = max(c(v1, v2)), SF, maxD = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{v1}{
%%     ~~Describe \code{v1} here~~
}
  \item{v2}{
%%     ~~Describe \code{v2} here~~
}
  \item{nbins}{
%%     ~~Describe \code{nbins} here~~
}
  \item{SF}{
%%     ~~Describe \code{SF} here~~
}
  \item{maxD}{
%%     ~~Describe \code{maxD} here~~
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
function (v1, v2, nbins = max(c(v1, v2)), SF, maxD = NULL) 
{
    if (is.null(maxD)) 
        maxD <- sum(abs(1:nbins - (nbins:1)))
    v1 <- as.factor(v1)
    v2 <- as.factor(v2)
    ints <- table(v1, v2)
    ridx <- as.numeric(rownames(ints))
    cidx <- as.numeric(colnames(ints))
    unions <- outer(table(v1), table(v2), "+") - ints
    unions[unions == 0] <- 1
    JSF <- SF[ridx, cidx] * (ints/unions)
    return(sum(JSF, na.rm = T)/maxD)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
