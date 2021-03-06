% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generateStream.R
\name{generate.static.stream}
\alias{generate.static.stream}
\title{Generate Static Data Stream}
\usage{
generate.static.stream(n = 1000, prop = 0.01,
  proptype = "proportional", stream.config = NULL, verbose = FALSE,
  method = "Rejection")
}
\arguments{
\item{n}{A vector containing \code{x} values, where the values corresponds
to the number of points for each step and \code{x} to the number of
steps.}

\item{prop}{Proportion of outliers in the hidden space.}

\item{proptype}{Type of the proportion of outliers. Value "proportional":
depend on the size of the empty space. Value "absolute": same absolute
proportion per subspace.}

\item{stream.config}{A stream configuration object. Should have been
generated with \code{nstep = 1}.}

\item{verbose}{Prints the number of the currently generated element if TRUE.}

\item{method}{Choose method of point generation. Can be "Rejection" or "Construction"}
}
\value{
An object of class stream, which is a List of 5 elements.
\itemize{
\item \code{data} contains the stream generated
\item \code{labels} contains the description of each point (\code{0} if the point
is not an outlier, or the subspace in which it is outlying as a string)
\item \code{n} the number of points at each step
\item \code{prop} the proportion of outliers in the hidden space
\item \code{stream.config} the associated stream configuration object (which is
valid only for static streams)
}
}
\description{
Generate a new synthetic multidimensional static data stream having the
desired properties.
}
\details{
The data is generated uniformly, except in certain subspaces where the data
is concentrated in particular dependencies (i.e. in the "Wall" dependency,
data concentrates on the axes, in a L-like shape). This should create spaces
with high dependency and space to observe hidden outliers. Note that the
proportion of outlier \code{prop} does not relate directly to the percentage
of outliers in the output stream. Since it corresponds to the probability of
a point, being ALREADY in the hidden space to stay where it is, the overall
proportion of outliers depends on the hidden space volume, which depends
on the number of subspaces and their margins. The greater the margin, the
bigger the hidden space.
}
\examples{
# Generate a stream with default parameters 
stream <- generate.static.stream()
# Generate a stream with custom configuration
stream.config <- generate.stream.config(dim=50, nstep=1) # nstep should be 1
stream <- generate.static.stream(n=1000, prop=0.05,
                                 stream.config=stream.config)
# Output stream results (to uncomment)
# output.stream(stream, "example")

}
\seealso{
\itemize{
\item \code{\link{generate.stream.config}}  : generate a stream.config file for a
dynamic or static stream
}
}
\author{
Edouard Fouché, \email{edouard.fouche@kit.edu}
}
