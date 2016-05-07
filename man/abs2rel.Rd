\name{abs2rel}
\alias{abs2rel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
abs2rel(qmix, idx)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{qmix}{
%%     ~~Describe \code{qmix} here~~
}
  \item{idx}{
%%     ~~Describe \code{idx} here~~
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
function (qmix, idx) 
{
    idx2d <- function(breaks, v) {
        idx <- matrix(NA, length(v), 2)
        for (i in 1:length(v)) {
            x <- sum(v[i] > breaks) + 1
            if (x > 1) 
                y <- v[i] - breaks[x - 1]
            else y <- v[i]
            idx[i, ] <- c(x, y)
        }
        return(idx)
    }
    breaks <- cumsum(sapply(qmix, ncol))
    idx <- idx2d(breaks, idx)
    return(idx)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
