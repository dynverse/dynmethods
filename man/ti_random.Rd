% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ti_container.R
\name{ti_random}
\alias{ti_random}
\title{Inferring a trajectory inference using Control: random}
\usage{
ti_random(dummy_param = 0.5)
}
\arguments{
\item{dummy_param}{Dummy parameter \cr
numeric; default: 0.5; possible values between 0 and 1}
}
\value{
The trajectory model
}
\description{
Will generate a trajectory using Control: random. This method was wrapped inside a \href{https://github.com/dynverse/dynmethods/tree/master/containers/random}{container}.
}