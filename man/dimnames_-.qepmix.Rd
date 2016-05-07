\name{dimnames<-.qepmix}
\alias{dimnames<-.qepmix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dimnames<-.qepmix(qm, value)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{qm}{
%%     ~~Describe \code{qm} here~~
}
  \item{value}{
%%     ~~Describe \code{value} here~~
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
function (qm, value) 
{
    len <- sum(sapply(qm, ncol))
    if (len != length(value)) 
        stop("qm number of samples and nms length differ")
    k <- 1
    for (i in 1:length(qm)) {
        n <- ncol(qm[[i]])
        colnames(qm[[i]]) <- value[k:(k + n - 1)]
        k <- k + n
    }
    return(qm)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
