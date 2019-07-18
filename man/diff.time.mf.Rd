\name{diff.time.mf}
\alias{diff.time.mf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
diff.time.mf(time1, time2, origin, units = c("auto", "secs", "mins", "hours", "days", "weeks"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{time1}{
%%     ~~Describe \code{time1} here~~
}
  \item{time2}{
%%     ~~Describe \code{time2} here~~
}
  \item{origin}{
%%     ~~Describe \code{origin} here~~
}
  \item{units}{
%%     ~~Describe \code{units} here~~
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
function (time1, time2, origin, units = c("auto", "secs", "mins", 
    "hours", "days", "weeks")) 
{
    if (missing(origin)) {
        time1 <- as.POSIXct(time1)
        time2 <- as.POSIXct(time2)
    }
    else {
        time1 <- as.POSIXct(time1, origin = origin)
        time2 <- as.POSIXct(time2, origin = origin)
    }
    z <- unclass(time1) - unclass(time2)
    attr(z, "tzone") <- NULL
    units <- match.arg(units)
    if (units == "auto") 
        units <- if (all(is.na(z))) 
            "secs"
        else {
            zz <- min(abs(z), na.rm = TRUE)
            if (!is.finite(zz) || zz < 60) 
                "secs"
            else if (zz < 3600) 
                "mins"
            else if (zz < 86400) 
                "hours"
            else "days"
        }
    switch(units, secs = .difftime(z, units = "secs"), mins = .difftime(z/60, 
        units = "mins"), hours = .difftime(z/3600, units = "hours"), 
        days = .difftime(z/86400, units = "days"), weeks = .difftime(z/(7 * 
            86400), units = "weeks"))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
