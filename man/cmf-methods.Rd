% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{cmf-methods}
\alias{cmf-methods}
\alias{plot.cmf}
\alias{screeplot.cmf}
\alias{summary.cmf}
\alias{print.cmf}
\title{Methods for cmf results lists}
\usage{
\method{plot}{cmf}(
  x,
  select,
  line = TRUE,
  labelSelected = TRUE,
  defaultColour = "#00008b",
  highlightColour = "#e2bd36",
  ...
)

\method{screeplot}{cmf}(x, topn, ...)

\method{summary}{cmf}(object, ...)

\method{print}{cmf}(x, ...)
}
\arguments{
\item{x}{cmf object}

\item{select}{optional selection vector of variables to show}

\item{line}{whether to show a line at the chosen cutoff}

\item{labelSelected}{whether to label the selected mediators in the plot, not
used when fewer than 20 bars are shown.}

\item{defaultColour}{the colour for the bars lower than the cutoff}

\item{highlightColour}{the colour for the bars of the selected mediators}

\item{...}{other arguments passed to barplot and summary}

\item{topn}{only show the top n mediators}

\item{object}{cmf object}
}
\description{
Plotting, summarising, and printing cmf objects
}
\examples{
# generate some data
dat <- generateMed(a = (1:10)/20, b = (1:10)/20)
res <- cmf(dat)
# screeplot of the result
screeplot(res)
# manhattan style plot the result
plot(res)


}
\seealso{
\code{\link{cmf}}
}
