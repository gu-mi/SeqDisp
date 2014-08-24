
##' @title Integrated likelihood and posteriors of epsilon and log(phi) for each of the m genes.
##' 
##' @details Based on il.prior.var.1, this function returns the integrated likelihood for all m genes (as a vector)
##' Also, it returns a vector of MAP estimator epsilon*, and a vector of MAP estimator log(phi).
##' 
##' @param prior.var  a number
##' @param prior.log.phi  a m-by-n matrix of log(phi0) pre-estimated
##' @param y read count matrix
##' @param s library size
##' @param x design matrix
##' @param beta0 regression coefficient vector
##' 
##' @return a list
##' 
##' @export
##' 
##' @author Yanming Di
##' 
il.prior.var = function(prior.var, prior.log.phi, y, s, x, beta0=rep(NA, dim(x)[2])) {
  
  m = dim(y)[1];
  
  l = numeric(m);
  epsilon.map = numeric(m);
  log.phi.map = numeric(m);
  
  for (i in 1:m){
    res = il.prior.var.1(prior.var, prior.log.phi[i,], y[i,], s, x, beta0);
    l[i] = res$l;
    epsilon.map = res$epsilon.map;
    log.phi.map = res$log.phi.map;
  }
  
  list(l=l, epsilon.map=epsilon.map, log.phi.map=log.phi.map);
}