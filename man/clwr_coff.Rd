% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clwr_coff.R
\name{clwr_coff}
\alias{clwr_coff}
\title{Eliminating the Unmeasured Confounders and Estimating
bidirectional causal effect and lagged causal effect}
\usage{
clwr_coff(
  data = data,
  y1 = "y1",
  y2 = "y2",
  y3 = "y3",
  x1 = "x1",
  x2 = "x2",
  x3 = "x3",
  centre = "centre"
)
}
\arguments{
\item{data}{an optional data frame containing the variables in the model.}

\item{y1}{the variable name of outcome at time t1}

\item{y2}{the variable name of outcome at time t2}

\item{y3}{the variable name of outcome at time t3}

\item{x1}{the variable name of exposure at time t1}

\item{x2}{the variable name of exposure at time t2}

\item{x3}{the variable name of exposure at time t3}

\item{centre}{the variable name of centre}
}
\value{
\code{est_xy} the casual effect estimation
and it's 95%CI of X at time t2 on Y at time t3

\code{est_xx} the casual effect estimation
and it's 95%CI of X at time t2 on X at time t3

\code{est_yx} the casual effect estimation
and it's 95%CI of Y at time t2 on X at time t3

\code{est_yy} the casual effect estimation
and it's 95%CI of Y at time t2 on Y at time t3
}
\description{
Eliminating the Unmeasured Confounders and Estimating
bidirectional causal effect and lagged causal effect
}
