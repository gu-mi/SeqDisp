
#' @title The Danio Rerio (zebrafish) RNA-Seq Dataset
#'
#' @description We use the zebrafish embryo model to study the innate immune response against Staphylococcus epidermidis. Therefore, we injected S. epidermidis into the yolk at 2 hpf and took samples at 5 days post injection. See References below.
#' 
#' @details Overview Design: This deep sequence study was designed to determine the gene expression profile by Staphylococcus epidermidis infection. RNA was isolated from embryos at 5 days post injection. Wildtypes zebrafish embryos were micro-injected into the yolk (2hpf) with 20 CFU of S. epidermdis O-47 mCherry bacteria suspended in PVP (Polyvinylpyrrolidone), or Non-injected as a control. After injections embryos were transferred into fresh egg water and incubated at 28 Celsius. At 5 days post injection 100-200 embryos per group were snap-frozen in liquid nitrogen, and total RNA was isolated using TRIZOL reagent. See References below.
#' 
#' @usage data(zebrafish)
#' @format A 27903 by 8 matrix of RNA-Seq read frequencies.
#' @name zebrafish
#'
#' @docType data
#'
#' @references See \url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42846} for more details.
#' 
#' @keywords datasets
#' 
#' @examples
#' 
#' library(SeqDisp)
#' data(zebrafish)
#' head(zebrafish)
#' dim(zebrafish)
#' 
#' # Control group only:
#' control = zebrafish[ ,seq(1,4)]
#' 
#' # Treatment group only:
#' treatment = zebrafish[ ,seq(5,8)]
#' 
NULL
