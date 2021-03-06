% Generated by roxygen2 (4.0.2): do not edit by hand
\name{sim.pwr.plot}
\alias{sim.pwr.plot}
\title{Preparing and plotting all dispersion models using ggplot2}
\usage{
sim.pwr.plot(sim.data.obj, type = "TPR-FDR", beta = 1, model.vec = NULL,
  title.note = "", legend = TRUE, legend.just = c(0, 1), N = 500,
  x.limit = c(0, 1), panel.label = "A")
}
\arguments{
\item{sim.data.obj}{output from sim.data.prepare().}

\item{type}{which type of plot to draw? Either "TPR-FDR", "N-TP" or "PREC"}

\item{beta}{a non-negative real value for calculating the general F_beta measure (default = 1: recall and precision are evenly weighted)}

\item{model.vec}{a character vector specifying the names of the NB dispersion models. Currently supported include
"NBP", "NBQ", "NBS", "Common", "Genewise", "Tagwise-Common", "Tagwise-Trend", "Trended", "QL", "QLShrink", "QLSpline" and "Genewise-HOA".}

\item{title.note}{extra data annotation in parentheses.}

\item{legend}{whether to put a legend on a single plot.}

\item{legend.just}{the legend.justification argument in ggplot2 controling the position of the legend.}

\item{N}{the number of top genes selected}

\item{x.limit}{specify the plotting range for the x-axis (TPR), for better visualization}

\item{panel.label}{the label on the lower-right corner of each plot (for journal figure legends).}
}
\value{
a ggplot object of the specified curves. Use print() to draw the specified curves.
}
\description{
Preparing and plotting all dispersion models using ggplot2
}
\details{
This function prepares meta-data and calls the following functions to obtain ggplot2 results:

   tpr.fdr.plot(): for TPR-FDR plot
   n.tp.plot(): for true positive plot (among top n genes)
   prec.plot(): for Precision-Recall plot
}
\author{
Gu Mi
}

