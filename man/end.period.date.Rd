\name{end.period.date}
\alias{end.period.date}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
end.period.date(date, interval)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{date}{
%%     ~~Describe \code{date} here~~
}
  \item{interval}{
%%     ~~Describe \code{interval} here~~
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
function (date, interval) 
{
    date.lt <- as.POSIXlt(date)
    switch(interval, y = {
        date.lt$mon = 11
        date.lt$mday = 31
        date = as.Date(date.lt)
    }, q = {
        date.lt$mon = (date.lt$mon\%/\%3 + 1) * 3\%\%12
        date.lt$mday = 1
        date.lt$year = date.lt$year + as.integer(date.lt$mon == 
            0)
        date = as.Date(date.lt) - 1
    }, m = {
        date.lt$mon = (date.lt$mon + 1)\%\%12
        date.lt$mday = 1
        date.lt$year = date.lt$year + as.integer(date.lt$mon == 
            0)
        date = as.Date(date.lt) - 1
    }, d = {
        date = as.Date(date)
    }, date = as.Date(date$lt))
    return(date)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
