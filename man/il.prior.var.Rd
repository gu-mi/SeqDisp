% Generated by roxygen2 (4.0.2): do not edit by hand
\name{il.prior.var}
\alias{il.prior.var}
\title{Integrated likelihood and posteriors of epsilon and log(phi) for each of the m genes.}
\usage{
il.prior.var(prior.var, prior.log.phi, y, s, x, beta0 = rep(NA, dim(x)[2]))
}
\arguments{
\item{prior.var}{a number}

\item{prior.log.phi}{a m-by-n matrix of log(phi0) pre-estimated}

\item{y}{read count matrix}

\item{s}{library size}

\item{x}{design matrix}

\item{beta0}{regression coefficient vector}
}
\value{
a list
}
\description{
Integrated likelihood and posteriors of epsilon and log(phi) for each of the m genes.
}
\details{
Based on il.prior.var.1, this function returns the integrated likelihood for all m genes (as a vector)
Also, it returns a vector of MAP estimator epsilon*, and a vector of MAP estimator log(phi).
}
\author{
Yanming Di
}

