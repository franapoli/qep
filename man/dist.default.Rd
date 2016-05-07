\name{dist.default}
\alias{dist.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dist.default(x, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{diag}{
%%     ~~Describe \code{diag} here~~
}
  \item{upper}{
%%     ~~Describe \code{upper} here~~
}
  \item{p}{
%%     ~~Describe \code{p} here~~
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
function (x, method = "euclidean", diag = FALSE, upper = FALSE, 
    p = 2) 
{
    if (!is.na(pmatch(method, "euclidian"))) 
        method <- "euclidean"
    METHODS <- c("euclidean", "maximum", "manhattan", "canberra", 
        "binary", "minkowski")
    method <- pmatch(method, METHODS)
    if (is.na(method)) 
        stop("invalid distance method")
    if (method == -1) 
        stop("ambiguous distance method")
    x <- as.matrix(x)
    N <- nrow(x)
    attrs <- if (method == 6L) 
        list(Size = N, Labels = dimnames(x)[[1L]], Diag = diag, 
            Upper = upper, method = METHODS[method], p = p, call = match.call(), 
            class = "dist")
    else list(Size = N, Labels = dimnames(x)[[1L]], Diag = diag, 
        Upper = upper, method = METHODS[method], call = match.call(), 
        class = "dist")
    .Call(C_Cdist, x, method, attrs, p)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
