################################################################################
## Copyright (C) 2014 Gu Mi <mig@stat.oregonstate.edu>
## 
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
## 02110-1301, USA
################################################################################

# Package Documentation
# 
# Author: Gu Mi
################################################################################

#' @title Evaluations of Negative Binomial Dispersion Methods for RNA Sequencing Data
#' 
#' @description An R package for evaluating the negative binomial dispersions methods commonly used in RNA-Seq data analysis
#' 
#' @details Details of the package
#' 
#' @name SeqDisp-package
#' @aliases SeqDisp-package SeqDisp
#' @docType package
#' 
#' @import NBPSeq edgeR QuasiSeq splines numDeriv plyr dplyr ggplot2 scales grid
#' 
#' @importFrom plyr llply
#' @importFrom plyr ldply
#' 
#' @author Yanming Di <diy@@stat.oregonstate.edu>, Gu Mi <neo.migu@@gmail.com>, Daniel W. Schafer
#' 
#' Maintainer: Gu Mi <https://github.com/gu-mi>
#'
#' @references See \url{https://github.com/gu-mi/SeqDisp/wiki/} for more details.
#' 
#' @keywords package 
#' 
NULL

.onLoad <- function(libname, pkgname){
  
  # suppress loading package messages
  suppressMessages(library(NBPSeq))
  suppressMessages(library(edgeR))
  suppressMessages(library(QuasiSeq))
  suppressMessages(library(splines))
  suppressMessages(library(numDeriv))
  suppressMessages(library(dplyr))
  suppressMessages(library(plyr))
  suppressMessages(library(ggplot2))
  suppressMessages(library(scales))
  suppressMessages(library(grid))
  #
  message("\n For updates of the SeqDisp package, please visit https://github.com/gu-mi/SeqDisp \n")
}  

