
##' @title Estimate extra variation in dispersion after a dispersion model has been fit
##' 
##' @description This function is similar to estimate.dispersion.var, except that it accepts log.phi directly
##' 
##' @param y read count matrix
##' @param s library size
##' @param x design matrix
##' @param log.phi pre-estimate dispersion matrix
##' @param lower  lower bound in \code{optimize}
##' @param upper  upper bound in \code{optimize}
##' @param hessian (logical) whether to calculate the Hessian matrix
##' 
##' @return an object
##' 
##' @export
##' 
##' @author Yanming Di
##' 
optimize.il.prior.var= function(y, s, x, log.phi, lower = 0, upper=200, hessian = FALSE) {
  
  ## if (!is.matrix(y)) y = matrix(y, 1, length(y));
  m = dim(y)[1];
  
  ll = function(prior.var) {
    l = 0;
    for (i in 1:m) {
      l = l + log(il.prior.var.1(prior.var, log.phi[i,], y[i,], s, x)$l);
      # il.prior.var.1: returns "l = il, epsilon.map = epsilon.map, log.phi.map = prior.log.phi + epsilon.map"
      # this loop returns just the likelihood "il"
      # "il" returned by il.prior.var.1 equals: exp(obj$objective) / sqrt(prior.var * d2l)
      # NOTE: d2l, the negative hessian, is calculated numerically --> better than analytically (if possible)?
    }
    l
  }
  
  # optimize the posterior log likelihood to get the sigma^2 that max the PLL (and hessian matrix if requested)
  # NOTE (mig): would it be better to have prior.var on the log scale, so its range stretches the real line (-inf, inf) ?
  # currently, prior.var is bounded by (0,inf)
  # if optimize on the log(var) scale, then exp(obj$maximum) may give a better estimate of sigma^2
  
  obj = optimize(ll, c(lower, upper), maximum=TRUE); 
  
  # if ask for hessian matrix, calculate numerically
  if (hessian) {
    obj$hessian = numDeriv::hessian(ll, obj$maximum);
  }
  
  obj
}