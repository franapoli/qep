\name{qep_scheme}
\alias{qep_scheme}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
qep_scheme(nbins, type = "uniform", data = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nbins}{
%%     ~~Describe \code{nbins} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (nbins, type = "uniform", data = NULL, ...) 
{
    pars <- list(...)
    if (type == "uniform") {
        if (!is.null(data)) 
            warning("data will be ignored by uniform scheme")
        scheme <- rep(nrow(data)/nbins, nbins)/nrow(data)
    }
    else if (type == "logistic") {
        if (!is.null(data)) {
            zeros <- apply(data, 2, function(x) sum(x > 0)/length(x))
            scheme <- lapply(zeros, function(x) qep_logistic(nbins, 
                pars$steep, pars$width, pars$baseline, x))
            scheme <- do.call(cbind, scheme)
        }
        else scheme <- qep_logistic(nbins, pars$steep, pars$width, 
            pars$baseline, pars$shift)
    }
    else if (type == "norm") {
        scheme <- dnorm(seq(-1, 1, length = nbins), sd = pars$sd)
        scheme <- scheme/sum(scheme)
    }
    else {
        stop("avaiable methds are: uniform, norm, logistic")
    }
    return(scheme)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
