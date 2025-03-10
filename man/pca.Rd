% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generic_functions.R
\name{pca}
\alias{pca}
\alias{pca.default}
\alias{pca.dfts}
\title{Generic pca of Data}
\usage{
pca(object, TVE = 1, ...)

\method{pca}{default}(object, ...)

\method{pca}{dfts}(object, TVE = 1, ...)
}
\arguments{
\item{object}{Data to compute principal component analysis on}

\item{TVE}{Numeric in [0,1] for the total variance explained. Can use this to
select only the required components.}

\item{...}{Additional information for \code{prcomp}.}
}
\value{
PCA data

Principal component data
}
\description{
Build on \code{prcomp}.
}
\examples{
pca(1:10)
pca(electricity)
}
