% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generateConfig.R
\name{replace.subspaces}
\alias{replace.subspaces}
\title{Replace a set of subspaces}
\usage{
replace.subspaces(dim, subspaces, indexes, allowOverlap = FALSE)
}
\arguments{
\item{dim}{Total number of dimensions from which to generate the subspaces.}

\item{subspaces}{A list of generated subspaces.}

\item{indexes}{Positions of the subspaces to replace.}

\item{allowOverlap}{States whether overlap between subspaces are allowed
(might be a bit more complex if TRUE). Note that full overlap and
subset relation are never allowed.}
}
\value{
The original set of subspaces, where the places at \code{indexes}
were modified
}
\description{
Resample compatible subspaces to replace a set of subspaces.
The new subspaces should not include or be included in another subspace from
the original set.
}
\examples{
subspaces <- list(c(1,2,3), c(3,4,5), c(7,8), c(11,12), c(15,16))
subspaces <- replace.subspaces(20, subspaces, indexes=c(2,4))

}
\author{
Edouard Fouché, \email{edouard.fouche@kit.edu}
}
