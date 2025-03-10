% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dfts.R
\name{minmax}
\alias{minmax}
\alias{max.dfts}
\alias{min.dfts}
\title{Max / Min for dfts Objects}
\usage{
\method{max}{dfts}(x, type = c("Obs", "intratime"), na.rm = TRUE, ...)

\method{min}{dfts}(x, type = c("Obs", "intratime"), na.rm = TRUE, ...)
}
\arguments{
\item{x}{A dfts object or data which can be automatically converted to that
format. See \code{\link[=dfts]{dfts()}}.}

\item{type}{String indicating if finding for observation ('Obs', default),
or for pointwise values ('intratime').}

\item{na.rm}{Boolean if NA values should be removed. Defaults to TRUE.}

\item{...}{Additional parameters to stats function}
}
\value{
A dfts object
}
\description{
Get the observation(s) with the min / max values. Selected as the observation
with the largest / smallest mean or the pointwise values.
}
\examples{
results <- max(electricity)
results <- min(electricity)
}
