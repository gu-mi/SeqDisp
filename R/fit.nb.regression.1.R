

##' Find mle of the dispersion and the regression coefficients in a NB regression model
##'
##' @title Find mle of the dispersion and the regression coefficients in a NB regression model
##' @param y 
##' @param s 
##' @param x 
##' @param beta0 
##' @param kappa 
##' @return a list
##' @author Yanming Di
##' 
fit.nb.regression.1 = function(y, s, x, beta0=rep(NA, dim(x)[2]), kappa=NA) {
  
  ## Find preliminary estiamtes of mu assuming phi=0.1. Will serve as
  ## initial values for the later irls algorithm.
  mustart = irls.nb.1(y, s, x, 0.1, beta0)$mu;
  
  if (is.na(kappa)) {
    
    ## Log likelihood of log(kappa)
    ll =  function(lkappa) {
      kappa = exp(lkappa);
      res = irls.nb.1(y, s, x, 1/kappa, beta0, mustart);
      sum(dnbinom(y, size =kappa, mu=res$mu, log=TRUE)); 
      ## log.likelihood.nb(kappa, res$mu, y);
    }
    
    res = optimize(ll, c(log(1e-20), log(1e20)), maximum=TRUE);
    kappa.hat = exp(res$maximum);
  } else {
    kappa.hat = kappa
  }
  
  phi.hat = 1/kappa.hat;
  
  res = irls.nb.1(y, s, x, phi.hat, beta0, mustart);
  beta.hat = res$beta;
  
  mu.hat = res$mu;
  ## l.hat = log.likelihood.nb(kappa.hat, mu.hat, y);
  l.hat = ll.nb(kappa.hat, mu.hat, y);
  
  v.hat = drop(mu.hat + phi.hat * mu.hat^2);
  n = length(y);
  n.pars = length(beta0) + 1;
  
  j.hat = matrix(NA, n.pars, n.pars);
  j.hat[-1,-1] = t(x) %*% diag(mu.hat^2 * (y/mu.hat^2 - (y+kappa.hat)/(mu.hat+kappa.hat)^2) - (y-mu.hat)*mu.hat/v.hat) %*% x;
  ## Be careful with the sign!
  j.hat[1,-1] = j.hat[-1,1] = - matrix((y - mu.hat) * mu.hat / (mu.hat+kappa.hat)^2, 1, n) %*% x;
  j.hat[1, 1] = -sum(trigamma(kappa.hat+y) - trigamma(kappa.hat) + 1/kappa.hat - 2/(kappa.hat+mu.hat) + (y+kappa.hat)/(mu.hat + kappa.hat)^2);
  
  ##  i.hat = t(x) %*% diag(mu.hat^2/v.hat) %*% x;
  ## l = log.likelihood.nb(kappa.hat, mu = mu.hat, y);
  l = ll.nb(kappa.hat, mu = mu.hat, y);
  
  ## See whether Poisson model fits better
  ## l.poisson  = sum(dpois(y, lambda = mu.hat, log=TRUE));
  ## eps = 1e-7;
  ## if (l.poisson + eps >= l & kappa.hat > 1e8) {
  ## l =  l.poisson;
  ##   kappa.hat = Inf;
  ##}
  ## l = sum(dnbinom(y, size = 1/phi, mu=res$mu, log=TRUE)); 
  
  list(kappa = kappa.hat, beta = beta.hat, mu = mu.hat, l=l.hat, j=j.hat);
}