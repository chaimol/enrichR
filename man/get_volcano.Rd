% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/volcano.R
\name{get_volcano}
\alias{get_volcano}
\title{get_volcano: get a volcano image}
\usage{
get_volcano(
  resdata,
  foldchange = 1,
  xrange = c(-10, 10),
  yrange = c(0, 100),
  legends = FALSE
)
}
\arguments{
\item{resdata}{a dataframe , fist col is Foldchange,second col is pvalue.}

\item{foldchange}{Default:1}

\item{xrange}{Default:c(-10:10)}

\item{yrange}{Default:c(0,100)}

\item{legends}{Ture or FALSE. Default:FALSE}
}
\value{
ggplot2 image
}
\description{
get_volcano: get a volcano image
}
\examples{
load("data/volcano.Rdata")
get_volcano(data1, foldchange = 1, xrange = c(-10, 10), yrange = c(0, 120), legends = TRUE)
get_volcano(data1, foldchange = 1.5)
get_volcano(data1)
}
