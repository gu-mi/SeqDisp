
#' @title The Homo Sapiens (human) RNA-Seq Dataset (Sequencing Depth = 5 Million)
#'
#' @description Motivation: RNA-seq is replacing microarrays as the primary tool for gene expression studies. Many RNA-seq studies have used insufficient biological replicates, resulting in low statistical power and inefficient use of sequencing resources. 
#' 
#' Results: We show the explicit trade-off between more biological replicates and deeper sequencing in increasing power to detect differentially expressed (DE) genes. In the human cell line MCF-7, adding more sequencing depth after 10M reads gives diminishing returns on power to detect DE genes, while adding biological replicates improves power significantly regardless of sequencing depth. We also propose a cost-effectiveness metric for guiding the design of large scale RNA-seq DE studies. Our analysis showed that sequencing less reads and perform more biological replication is an effective strategy to increase power and accuracy in large scale differential expression RNA-seq studies, and provided new insights into efficient experiment design of RNA-seq studies. See References below.
#' 
#' @details Overall Design: Treatment (10nM E2 treatment for 24h) and control MCF7 cells are both replicated 7 times, and collected for mRNA-seq. Reads are then subsampled for statistical analysis. See References below.
#' 
#' @usage data(human5)
#' @format A 22336 by 14 matrix of RNA-Seq read frequencies.
#' @name human5
#'
#' @docType data
#'
#' @references Liu Y, Zhou J, White KP. RNA-seq differential expression studies: more sequence or more replication? 
#' Bioinformatics 2014 Feb 1;30(3): 301-4. PMID: 24319002
#' 
#' See \url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE51403} for more details.
#' 
#' @keywords datasets
#' @examples
#' 
#' library(SeqDisp)
#' 
#' data(human5)
#' head(human5)
#' dim(human5)
#' 
#' # Control group only:
#' control = human5[ ,seq(1,7)]
#' 
#' # Treatment group only:
#' treatment = human5[ ,seq(8,14)]
#' 
NULL