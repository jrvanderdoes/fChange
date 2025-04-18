% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pca_examination.R
\name{projection_model}
\alias{projection_model}
\title{Projection-based functional data model}
\usage{
projection_model(
  X,
  TVE = 0.95,
  model = c("ets", "arima"),
  n.ahead = 0,
  alpha = 0.05,
  check.cp = TRUE,
  frequency = 1,
  ...
)
}
\arguments{
\item{X}{A dfts object or data which can be automatically converted to that
format. See \code{\link[=dfts]{dfts()}}.}

\item{TVE}{Numeric in [0,1] for the total variance explained to select
number of PCA components to use to model the data.}

\item{model}{String to indicate method to model components, either
"ets" or "arima".}

\item{n.ahead}{Number of observations to forecast.}

\item{alpha}{Significance in [0,1] for intervals when forecasting.}

\item{check.cp}{Boolean which indicates if the errors should be checked for
change points to change forecasts and plots.}

\item{frequency}{Numeric for seasonal frequency when component is made a ts
object for the models.}

\item{...}{Additional information to pass into pca, change (if
\code{check.cp=TRUE}), and plot.}
}
\value{
List with the following elements:
\itemize{
\item fit: dfts object for fit.
\item forecast_plot: plot of the data with any forecasted values.
\item residuals: dfts object for residuals from the fit.
\item changes: vector of any changes when using \code{detect.cp}.
\item component_models: modeled PCs from the data.
\item component_true: true data constucted via the PCs.
\item parameters: list with fit parameters like pcs, TVE, model, and n.ahead.
}
}
\description{
Model and forecast functional data using a Hyndman and Ullah projection-based
model.
}
\examples{
result <- projection_model(dfts(electricity$data[,50:100]),
 n.ahead=1, TVE=0.1, check.cp=FALSE)
}
\references{
Hyndman, R. J., & Shahid Ullah, M. (2007). Robust forecasting
of mortality and fertility rates: A functional data approach. Computational
Statistics & Data Analysis, 51(10), 4942-4956.
https://doi.org/10.1016/j.csda.2006.07.028
}
