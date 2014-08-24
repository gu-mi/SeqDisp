
#' @title The Arabidopsis RNA-Seq Dataset
#'
#' @description The goal of this study is to identify the targets of RVE8, a MYB-like transcription factor involved in the circadian clock in Arabidopsis. Analysis of 7 days old rve8-1 RVE8::RVE8:GR and rve8-1 seedlings treated with dexamethasone or mock identified genes responsive to RVE8 induction.The RVE8 up-regulated genes are enriched for evening-phased genes while the down-regulated genes are enriched for a morning phase. This study reveals that RVE8 is a master regulator of circadian gene expression in Arabidopsis. See References below.
#' 
#' @details Overview Design: Transgenic line rve8-1 RVE8::RVE8:GR and rve8-1 treated with DEX or mock with three biological replicates each, 12 samples in total. See References below.
#' 
#' @usage data(arabidopsis)
#' @format A 21567 by 12 matrix of RNA-Seq read frequencies.
#' @name arabidopsis
#'
#' @docType data
#'
#' @references See \url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE38879} for more details.
#' 
#' @keywords datasets
#' 
#' library(SeqDisp)
#' data(arabidopsis)
#' head(arabidopsis)
#' dim(arabidopsis)
#' 
#' # RVE8:GR_mock group only:
#' RVE8_GR_mock = arabidopsis[ ,1:3]
#' 
#' # RVE8:GR_DEX group only:
#' RVE8_GR_DEX = arabidopsis[ ,4:6]
#'
#' # rve8_mock group only:
#' rve8_mock = arabidopsis[ ,7:9]
#'  
#' # rve8_DEX group only:
#' rve8_DEX = arabidopsis[ ,10:12]
#' 
#' 
NULL
