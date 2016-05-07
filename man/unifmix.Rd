\name{unifmix}
\alias{unifmix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
unifmix(qm, v, subslice = F, na.rm = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{qm}{
%%     ~~Describe \code{qm} here~~
}
  \item{v}{
%%     ~~Describe \code{v} here~~
}
  \item{subslice}{
%%     ~~Describe \code{subslice} here~~
}
  \item{na.rm}{
%%     ~~Describe \code{na.rm} here~~
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
function (qm, v, subslice = F, na.rm = T) 
{
    if (!subslice) {
        v2 <- sapply(v, function(x) rel2abs(qmix, x, 1:ncol(qmix[[x]])))
        if (length(v) > 1) 
            v <- do.call(c, v2)
        else v <- v2
    }
    idx <- abs2rel(qm, v)
    qeps <- unique(idx[, 1])
    k <- 1
    if (na.rm) {
        nonNAset <- rownames(qmix[[qeps[1]]])
        if (length(qeps) > 1) 
            for (i in 2:length(qeps)) nonNAset <- intersect(nonNAset, 
                rownames(qm[[qeps[i]]]))
    }
    else nonNAset <- attr(qm, "genes")
    resmat <- list()
    for (i in 1:length(qeps)) {
        q <- qeps[i]
        if (na.rm) {
            sub <- qm[[q]][nonNAset, idx[idx[, 1] == q, 2], drop = F]
        }
        else sub <- qm[[q]][, idx[idx[, 1] == q, 2], drop = F]
        resmat[[i]] <- qep(sub)
    }
    return(qepmix(resmat))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
