\name{lag.num}
\alias{lag.num}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
lag.num(x.lag, period, unit)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x.lag}{
%%     ~~Describe \code{x.lag} here~~
}
  \item{period}{
%%     ~~Describe \code{period} here~~
}
  \item{unit}{
%%     ~~Describe \code{unit} here~~
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
function (x.lag, period, unit) 
{
    if (is.numeric(x.lag) == T && is.atomic(x.lag) == T && length(x.lag) == 
        1L) {
        return(x.lag)
    }
    multiplier <- as.double(substr(x.lag, start = 1, stop = nchar(x.lag) - 
        1))
    if (is.na(multiplier)) {
        stop("The description of lags cannot be recognized. The format should be 3m, 1q, etc")
    }
    ndaysPerYear <- 264
    ndaysPerQuarter <- 66
    ndaysPerMonth <- 22
    nhoursPerDay <- 8
    if (substr(x.lag, start = nchar(x.lag), stop = nchar(x.lag)) == 
        "y") {
        multiplier <- multiplier * ndaysPerYear
    }
    if (substr(x.lag, start = nchar(x.lag), stop = nchar(x.lag)) == 
        "q") {
        multiplier <- multiplier * ndaysPerQuarter
    }
    if (substr(x.lag, start = nchar(x.lag), stop = nchar(x.lag)) == 
        "m") {
        multiplier <- multiplier * ndaysPerMonth
    }
    if (substr(x.lag, start = nchar(x.lag), stop = nchar(x.lag)) == 
        "d") {
        multiplier <- multiplier * 1
    }
    if (substr(x.lag, start = nchar(x.lag), stop = nchar(x.lag)) == 
        "h") {
        multiplier <- multiplier/nhoursPerDay
    }
    if (substr(x.lag, start = nchar(x.lag), stop = nchar(x.lag)) == 
        "s") {
        multiplier <- multiplier/(nhoursPerDay * 60 * 60)
    }
    if (unit == 1) {
        x.lag <- round(multiplier/(ndaysPerYear * period))
    }
    if (unit == 2) {
        x.lag <- round(multiplier/(ndaysPerMonth * period))
    }
    if (unit == 3) {
        x.lag <- round(multiplier/period)
    }
    if (unit == 4) {
        x.lag <- round(multiplier/(period/nhoursPerDay))
    }
    if (unit == 5) {
        x.lag <- round(multiplier/(period/nhoursPerDay/60))
    }
    if (unit == 6) {
        x.lag <- round(multiplier/(period/nhoursPerDay/60/60))
    }
    return(x.lag)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
