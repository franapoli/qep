\name{qep_checknames}
\alias{qep_checknames}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
qep_checknames(qdata)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{qdata}{
%%     ~~Describe \code{qdata} here~~
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
function (qdata) 
{
    gnames <- rownames(qdata)
    if (class(qdata) != "matrix") 
        stop("qdata must be an object of class matrix.")
    if (is.null(gnames)) 
        stop(paste0("Please provide gene names as rownames ", 
            "of qdata or through the gnames parameter."))
    if (any(duplicated(gnames))) 
        stop("Gene names must be unique.")
    if (any(is.na(gnames))) 
        stop("Gene names can not be NA")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
