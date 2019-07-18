\name{data.freq}
\alias{data.freq}
\title{
data.freq
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
data.freq(DateVec)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DateVec}{
%%     ~~Describe \code{DateVec} here~~
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
function (DateVec) 
{
    DateDiff <- as.matrix(diff(DateVec))
    modeUse = mode(DateDiff[, 1])$dataMode
    if (modeUse >= 1) {
        period <- modeUse
        unit <- 1
        return(list(period = period, unit = unit))
    }
    modeUse <- mode(DateDiff[, 2])$dataMode
    mask <- isTRUE(modeUse < 0)
    modeUse[mask] <- modeUse[mask] + 12
    if (modeUse >= 1) {
        period <- modeUse
        unit <- 2
        return(list(period = period, unit = unit))
    }
    modeUse <- mode(DateDiff[, 3])$dataMode
    mask <- isTRUE(modeUse < 0)
    modeUse[mask] <- modeUse[mask] + 30
    if (modeUse >= 1) {
        period <- modeUse
        unit <- 3
        return(list(period = period, unit = unit))
    }
    modeUse <- mode(DateDiff[, 4])$dataMode
    mask <- isTRUE(modeUse < 0)
    modeUse[mask] <- modeUse[mask] + 24
    if (modeUse >= 1) {
        period <- modeUse
        unit <- 4
        return(list(period = period, unit = unit))
    }
    modeUse <- mode(DateDiff[, 5])$dataMode
    mask <- isTRUE(modeUse < 0)
    modeUse[mask] <- modeUse[mask] + 60
    if (modeUse >= 1) {
        period <- modeUse
        unit <- 5
        return(list(period = period, unit = unit))
    }
    elapse <- diff.time.mf(DateVec[2:dim(DateVec)[1], 6], DateVec[1:dim(DateVec)[1] - 
        1, 6], origin = "1970-01-01", units = "secs")
    period <- mean(elapse)
    unit <- 6
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
