\name{is.seqelist}
\alias{is.eseq}
\alias{is.seqe}
\alias{is.seqelist}
\alias{print.eseq}
\alias{print.seqelist}
\alias{plot.eseq}
\alias{plot.seqelist}
\alias{levels.seqelist}
\alias{levels.eseq}
\alias{[.seqelist}
\alias{str.eseq}
\alias{str.seqelist}
\alias{summary.seqelist}
\alias{summary.eseq}
\alias{print.eseqsum}
%\alias{Math.eseq}
%\alias{Math.seqelist}
%\alias{Ops.eseq}
%\alias{Ops.seqelist}
%\alias{Summary.eseq}
%\alias{Summary.seqelist}


\title{Event sequence object}
\description{
  TraMineR uses an internal event sequence object for all its operations
  with event sequences such as mining frequent subsequences. 
  
  Function \code{is.seqelist} checks wether the argument is an event
  sequence object while \code{is.eseq} checks wether the argument is a
  single element of an event sequence object.
  %\code{as.seqelist} transform a list of event sequences into an internal event sequence object.

  is.seqe() has been removed in TraMineR 2.0-4, use is.eseq() instead.
  
  There are subsetting, \code{print}, \code{plot}, \code{levels}, and \code{str} methods that apply to event sequence objects and their elements. The \code{summary} method applies to \code{seqelist} objects and returns number of sequences, size of the alphabet, and summary values of number of transitions per sequence and sequence durations.

}
\usage{
is.eseq(eseq, s)
is.seqelist(eseq, s)
\method{print}{eseq}(x, quote=FALSE, ...)
\method{print}{seqelist}(x, quote=FALSE, ...)
\method{plot}{eseq}(x, type="pc", ...)
\method{plot}{seqelist}(x, type="pc", ...)
\method{levels}{eseq}(x, ...)
\method{levels}{seqelist}(x, ...)
\method{str}{eseq}(object, ...)
\method{str}{seqelist}(object, ...)
\method{summary}{eseq}(object, ...)
\method{summary}{seqelist}(object, ...)
}
\arguments{
  \item{eseq}{Sequence object}
  \item{x, object}{\code{eseq} or \code{seqelist} object}
  \item{quote}{logical, should strings be printed with surrounding quotes?}
  \item{type}{string, currently only \code{"pc"} (parallel coordinate plot) is supported}
  \item{...}{further arguments to other functions or methods.}
  \item{s}{Deprecated. Use \code{eseq} instead.}
}
\value{The \code{summary} method applies to \code{seqelist} objects and returns number of sequences, size of the alphabet, and summary values of sequence durations (\code{dur}) and number of transitions (\code{ntrans}) per sequence. It also returns transparently the vectors \code{$dur} of sequence durations and \code{$ntrans} of number of transitions. When applied to an element of a \code{seqelist}, it returns the values of \code{ntrans} and \code{dur} of the selected element and \code{nalph}, the size of the alphabet.
}
\seealso{\code{\link{seqpcplot}}, \code{\link{as.character.seqelist}}, and \code{\link{as.character.eseq}}, the latter two for character string representations of event sequence objects}
\examples{
data(ex1)
sqex1 <- seqdef(ex1[,1:13])
esqex1 <- seqecreate(sqex1)
## esqex1 is a list of event sequences
is.seqelist(esqex1)
## alphabet of the event sequence object
levels(esqex1)
## can also be extracted with
alphabet(esqex1)

is.eseq(esqex1[[1]])

(res <- summary(esqex1))
res$ntrans
res$dur

summary(esqex1[[3]])
}

\author{Matthias Studer and Gilbert Ritschard}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Event sequences}
\keyword{internal}
