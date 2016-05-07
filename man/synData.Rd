\name{synData}
\alias{synData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
synData(nobsuniv, nobs, nconds, nbins, meansd = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nobsuniv}{
%%     ~~Describe \code{nobsuniv} here~~
}
  \item{nobs}{
%%     ~~Describe \code{nobs} here~~
}
  \item{nconds}{
%%     ~~Describe \code{nconds} here~~
}
  \item{nbins}{
%%     ~~Describe \code{nbins} here~~
}
  \item{meansd}{
%%     ~~Describe \code{meansd} here~~
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
function (nobsuniv, nobs, nconds, nbins, meansd = 0.05) 
{
    genome <- paste0("O", 1:nobsuniv)
    rawdata <- list()
    qdata <- list()
    for (i in 1:length(nobs)) {
        nr <- nobs[[i]]
        nc <- nconds[[i]]
        m <- matrix(rnorm(nr * nc, 0.5), nr, nc)
        m <- apply(m, 2, function(x) x + rnorm(1, sd = meansd))
        rawdata[[i]] <- m
        rownames(rawdata[[i]]) <- sample(genome, nr)
        colnames(rawdata[[i]]) <- paste0("P", i, "_C", 1:nc)
        s <- logistic.ep(nbins, steep = 0.5, width = 0, baseline = 0.1, 
            zero = rawdata[[i]])
        qdata[[i]] <- quantize.ep(rawdata[[i]], s)
    }
    qmix <- qepmix(qdata)
    return(list(rawdata = rawdata, qep = qdata))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
