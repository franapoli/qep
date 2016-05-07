\name{qep}
\alias{qep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
qep(qdata, gnames = rownames(qdata), cndnames = colnames(qdata))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{qdata}{
%%     ~~Describe \code{qdata} here~~
}
  \item{gnames}{
%%     ~~Describe \code{gnames} here~~
}
  \item{cndnames}{
%%     ~~Describe \code{cndnames} here~~
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
function (qdata, gnames = rownames(qdata), cndnames = colnames(qdata)) 
{
    if (!is.integer(qdata)) 
        stop("qdata must be of class integer.")
    if (any(is.na(qdata))) 
        stop("qdata can not contain NAs")
    qep_checknames(qdata)
    dimnames(qdata) <- list(genes = gnames, conditions = cndnames)
    return(structure(qdata, class = "qep"))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
