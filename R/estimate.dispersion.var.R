
##' @title Estimate extra variation in dispersion after a dispersion model has been fit (for NBPSeq)
##' 
##' @details The fitted dispersion models may come from the NBPSeq packages. This function are used for both simulation
##' (if we refit the dispersion model; see sim.dr) and real RNA-Seq data studies.
##'
##' @param nb.data output from \code{link{prepare.nb.data}}
##' @param dispersion output from \code{estimate.dispersion} (NBPSeq)
##' @param x model matrix, should be the same one used to estimate \code{dipserison}
##' @param subset.rows  the row indices for subsetting the read count matrix
##' @param lower  lower bound in \code{optimize}
##' @param upper  upper bound in \code{optimize}
##' @param hessian (logical) whether to calculate the Hessian matrix
##'
##' @return an object
##' 
##' @export
##' 
##' @author Yanming Di, Gu Mi
##' 
estimate.dispersion.var = function(nb.data, dispersion, x, subset.rows = 1:m, lower = 0, upper = 200, hessian = FALSE) {
  
  y = nb.data$counts
  s = nb.data$eff.lib.sizes
  m = dim(y)[1]
  row.ids = (1:m)[subset.rows]
  
  # since NBPSeq and edgeR report dispersion estimates differently, we need to be careful in extracting the dispersions!
  # NBPSeq:
  log.phi = log(dispersion$estimates)  # this log.phi is the theta0, pre-estimated by some model (e.g. NBQ)
  
  # posterior log likelihood of prior.var
  ll = function(prior.var) {
    l = 0
    for (i in row.ids) {
      l = l + log(il.prior.var.1(prior.var, log.phi[i,], y[i,], s, x)$l)  # adds on the log likelihood for each gene
      # il.prior.var.1: returns "l = il, epsilon.map = epsilon.map, log.phi.map = prior.log.phi + epsilon.map"
      # this loop returns just the likelihood "il"
      # "il" returned by il.prior.var.1 equals: exp(obj$objective) / sqrt(prior.var * d2l)
      # NOTE: d2l, the negative hessian, is calculated numerically
    }
    l
  }
  
  # optimize the posterior log likelihood to get the sigma^2 that max the PLL (and hessian matrix if requested)
  obj = optimize(ll, c(lower, upper), maximum=TRUE) 
  
  # if ask for hessian matrix, calculate numerically
  if (hessian) {
    obj$hessian = numDeriv::hessian(ll, obj$maximum)
  }
  
  obj
}


##' @title Estimate extra variation in dispersion after a dispersion model has been fit (for edgeR)
##' 
##' @details The fitted dispersion models may come from the edgeR packages. This function are used for both simulation
##' (see sim.dr) and real RNA-Seq data studies.
##'
##' @param nb.data output from \code{link{prepare.nb.data}}
##' @param dispersion dispersion matrix of dimension m-by-n (edgeR)
##' @param x model matrix, should be the same one used to estimate \code{dipserison}
##' @param subset.rows  the row indices for subsetting the read count matrix
##' @param lower  lower bound in \code{optimize}
##' @param upper  upper bound in \code{optimize}
##' @param hessian (logical) whether to calculate the Hessian matrix
##'
##' @return an object
##' 
##' @export
##' 
##' @author Yanming Di, Gu Mi
##' 
estimate.dispersion.var.edgeR = function(nb.data, dispersion, x, subset.rows = 1:m, lower = 0, upper = 200, hessian = FALSE) {
  
  y = nb.data$counts
  s = nb.data$eff.lib.sizes
  m = dim(y)[1]
  row.ids = (1:m)[subset.rows]
  
  # since NBPSeq and edgeR report dispersion estimates differently, we need to be careful in extracting the dispersions!
  # edgeR:
  log.phi = log(dispersion)
  
  # posterior log likelihood of prior.var
  ll = function(prior.var) {
    l = 0
    for (i in row.ids) {
      l = l + log(il.prior.var.1(prior.var, log.phi[i,], y[i,], s, x)$l)  # adds on the log likelihood for each gene
      # il.prior.var.1: returns "l = il, epsilon.map = epsilon.map, log.phi.map = prior.log.phi + epsilon.map"
      # this loop returns just the likelihood "il"
      # "il" returned by il.prior.var.1 equals: exp(obj$objective) / sqrt(prior.var * d2l)
      # NOTE: d2l, the negative hessian, is calculated numerically
    }
    l
  }
  
  # optimize the posterior log likelihood to get the sigma^2 that max the PLL (and hessian matrix if requested)
  obj = optimize(ll, c(lower, upper), maximum=TRUE) 
  
  # if ask for hessian matrix, calculate numerically
  if (hessian) {
    obj$hessian = numDeriv::hessian(ll, obj$maximum)
  }
  
  obj
}