\name{as.character.eseq}
%\alias{str.eseq}
%\alias{str.seqelist}
\alias{as.character.eseq}
\alias{as.character.seqelist}

\title{String representation of an event sequence object}
\description{
  Generic \code{as.character} methods to get a string representation of an event sequence \code{seqelist} object and  of an element \code{eseq} of an event sequence object.
  %These are the generic \code{str} and methods for objects of class (single event sequence) and \code{seqelist} (list of event sequences).
}
\usage{
\method{as.character}{seqelist}(x, ...)
\method{as.character}{eseq}(x, ...)
}

\arguments{
  \item{x}{Event sequence \code{seqelist} or \code{esqe} object to be coerced}
  \item{...}{Further arguments passed to other methods}
}

\details{
   String representations have the following form:
  \code{time-(e1,e2,...)-elapsedtime-(e3,...)...}
  where \code{time} is the time elapsed from start to the first set of simultaneous events \code{(e1,e2,...)} and \code{elapsedtime} the time between two consecutive sets of simultaneous events. Time is only displayed when events are time stamped.
  %Use \code{str} (\code{str.seqelist}) to get a single string for a whole list, and \code{as.charcater} (\code{as.character.seqelist}) to get a vector of strings.
  %\code{str.eseq} and \code{as.character.eseq} have mainly an internal purpose (invoked by print methods for example).
}

\seealso{\code{\link{seqecreate}} and \code{\link{print.seqelist}} for other methods applying to event sequence objects.}
\examples{
data(actcal.tse)
actcal.eseq <- seqecreate(actcal.tse)
head(as.character(actcal.eseq))
}

\author{Matthias Studer (with Gilbert Ritschard for the help page)}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

\keyword{internal}
\keyword{Event Sequences}
