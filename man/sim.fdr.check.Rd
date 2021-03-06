% Generated by roxygen2 (4.0.2): do not edit by hand
\name{sim.fdr.check}
\alias{sim.fdr.check}
\title{Simulate datasets to check FDR controls for different NB dispersion models}
\usage{
sim.fdr.check(model = "Genewise", m = 5000, n = 6, mean = 6.5,
  sd = 1.5, perc.de = 0.05, FDR = 0.1, seed.vec = NULL)
}
\arguments{
\item{model}{NB dispersion model name}

\item{m}{number of genes}

\item{n}{number of samples}

\item{mean}{the mean of a normal distribution}

\item{sd}{the standard deviation of a normal distribution}

\item{perc.de}{percentage of DE genes}

\item{FDR}{nominal FDR threshold specified below which we call a DE gene}

\item{seed.vec}{random number generator seed for reproducibility}
}
\value{
a data frame with variables to be used for plotting
}
\description{
We simulate read counts of two groups. The (log) mean is determined by normal distribution,
and the percentage of DE genes can be specified, with half up- and half down-regulated. Nominal FDR can be
specified, and we calculate the (empirical) true FDR and then compare.
}
\examples{
library(SeqDisp)
library(dplyr)

# check three NB dispersion models with two levels of percentage of DE
a1 = sim.fdr.check(model="Genewise", m=5000, n=6, mean=6.5, sd=1.5, perc.de=0.2, FDR=0.1, seed.vec=c(1,2,3))
b1 = sim.fdr.check(model="Trended", m=5000, n=6, mean=6.5, sd=1.5, perc.de=0.2, FDR=0.1, seed.vec=c(1,2,3))
c1 = sim.fdr.check(model="Tagwise-Trend", m=5000, n=6, mean=6.5, sd=1.5, perc.de=0.2, FDR=0.1, seed.vec=c(1,2,3))
a2 = sim.fdr.check(model="Genewise", m=5000, n=6, mean=6.5, sd=1.5, perc.de=0.05, FDR=0.1, seed.vec=c(1,2,3))
b2 = sim.fdr.check(model="Trended", m=5000, n=6, mean=6.5, sd=1.5, perc.de=0.05, FDR=0.1, seed.vec=c(1,2,3))
c2 = sim.fdr.check(model="Tagwise-Trend", m=5000, n=6, mean=6.5, sd=1.5, perc.de=0.05, FDR=0.1, seed.vec=c(1,2,3))

# combine the results
res = rbind(a1,b1,c1,a2,b2,c2)
tb_df = tbl_df(data=res)
tb_gb = group_by(x=tb_df, model)
tb_final = mutate(tb_gb, ratio = actual.fdr/nominal.fdr)

# graphical comparison using ggplot2
p = ggplot(data = tb_final, aes(x=model, y=ratio, colour=perc.de, shape=perc.de)) +
geom_boxplot(aes(colour=perc.de, shape=perc.de), outlier.size = 0.5, size=0.2) +
geom_hline(yintercept=1, colour="Black", size=0.2, linetype=2) +
theme_bw() +
scale_x_discrete("") +
scale_y_continuous("Actual FDR / Nominal FDR", limits = c(0, max(tb_final[ ,"ratio"]))) +
scale_colour_discrete(name =  "\% DE genes") +
scale_shape_discrete(name = "\% DE genes") +
theme(plot.title = element_text(face="bold", size=6),
# axis labels
axis.title.x = element_text(face="bold", size=5),
axis.title.y = element_text(face="bold", size=5, angle=90),
# axis tick labels
axis.text.x = element_text(angle=0, vjust=0.5, size=5),
axis.text.y = element_text(size=5),
axis.ticks.length = unit(0.5, "mm"),
legend.justification=c(1,1),
legend.position=c(1,1),
legend.title = element_text(size=6),
legend.key = element_blank(),
legend.text = element_text(size=4)
)
print(p)
}
\author{
Gu Mi
}

