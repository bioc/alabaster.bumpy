% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/readBumpyDataFrameMatrix.R
\name{readBumpyDataFrameMatrix}
\alias{readBumpyDataFrameMatrix}
\alias{loadBumpyDataFrameMatrix}
\title{Read a BumpyDataFrameMatrix from disk}
\usage{
readBumpyDataFrameMatrix(path, metadata, ...)
}
\arguments{
\item{path}{String containing a path to a directory, itself created using the \code{\link{saveObject}} method for \linkS4class{BumpyDataFrameMatrix} objects.}

\item{metadata}{Named list of metadata for this object, see \code{\link{readObjectFile}} for details.}

\item{...}{Further arguments passed to internal \code{\link{altReadObject}} calls.}
}
\value{
A \linkS4class{BumpyDataFrameMatrix} object.
}
\description{
Read a \linkS4class{BumpyDataFrameMatrix} from its on-disk representation.
}
\examples{
# Mocking up a BumpyMatrix.
library(BumpyMatrix)
library(S4Vectors)
df <- DataFrame(x=runif(100), y=runif(100))
f <- factor(sample(letters[1:20], nrow(df), replace=TRUE), letters[1:20])
out <- S4Vectors::split(df, f)
mat <- BumpyMatrix(out, c(5, 4))

# Saving it:
tmp <- tempfile()
saveObject(mat, tmp)

# Reading it:
readBumpyDataFrameMatrix(tmp)

}
\author{
Aaron Lun
}
