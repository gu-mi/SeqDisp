

##' This function finds the MAP estimator of log(phi) and computes
##' the integrated likelihood (IL) of the variance of the prior
##' distribution on log(phi) by Laplace approximation.   
##'
##' Note that the components of log(phi.prior) can be different, but
##' it is assumed that
##'
##'   epsilon = log(phi) - prior.log.phi
##'
##' is the same for all components and epsilon has a (prior) distribution
##'
##'   N(0, prior.var).
##'
##' Given phi, the components of y follow a NB distribution with
##' means mu = s exp(x' beta) and dispersion phi.
##'
##' The integrated likelihood of \code{prior.var} is defined as the integration of
##'
##'     Lp(epsilon; y) Prior(epsilon)
##'
##' over epsilon, where Lp(epsilon; y) is the profile likelihood of epsilon.
##'
##' For given epsilon, we can determine the value of phi=phi(epsilon),
##' and find the MLE of beta. The profile likelihood of epsilon is the
##' likelihood of (epsilon, beta) computed at epsilon and the MLE of
##' beta given epsilon: l(phi(epsilon), beta_MLE(phi)).
##'
##' @title Compute integrated likelihood of the variance of the
##' prior distribution of log dispersion
##'
##' @param prior.var a number, the variance of the prior distribution of log(phi) --> sigma^2
##' @param prior.log.phi a number or an n-vector, the mean(s) of the prior distribution of log(phi)  --> theta0
##' @param y a n-vector of NB counts.
##' @param s a n-vector of library sizes.
##' @param x a n by p design matrix.
##' @param beta0 a p-vector specifying the known and unknown
##' components of beta, the regression coefficients. NA values
##' indicate unknown components and non-NA values specify the values
##' of the known components. The default is that all components of
##' beta are unknown.
##' 
##' @return a list containing the following components:
##'   \item{l}{the integrated likelihood of \code{prior.var}}
##'   \item{epsilon.map}{a scalar, MAP estimate of epsilon = log.phi - prior.log.phi}
##'   \item{log.phi.map}{a scalar or a vector, MAP estimate of log.phi}
##'   
##' @author Yanming Di
##' 
##' @export
##'  
il.prior.var.1 = function(prior.var, prior.log.phi, y, s, x, beta0=rep(NA, dim(x)[2])) {
  
  ## Find preliminary estimates of mu assuming phi=0.1. Will serve as
  ## initial values for the later irls algorithm.
  #args(irls.nb.1)
  mustart = irls.nb.1(y, s, x, phi=0.1, beta0)$mu;
  
  ## Posterior log likelihood of epsilon
  ll =  function(epsilon) {
    kappa = exp(-(prior.log.phi + epsilon));  # exp(-theta) = exp(-log(phi)) = 1/phi --> phi = 1/kappa (used next)
    res = irls.nb.1(y, s, x, phi=1/kappa, beta0, mustart);
    l = sum(dnbinom(y, size = kappa, mu=res$mu, log=TRUE));   # this is the l_p(epsilon;y_i), log profile likelihood of epsilon
    l - epsilon^2/ (2 * prior.var)
    ## log.likelihood.nb(kappa, res$mu, y);
    
    ## NOTE (mig): this is the g(epsilon;y_i); next we need to find epsilon* that maximize it (epsilon* = epsilon.map)
  }
  
  ## Find the MAP estimate of epsilon = log.phi - prior.log.phi
  
  ## FIXME: these bounds might be too wide (for range of epsilon ~ N(0, sigma^2))
  lower = log(1e-10);
  upper = log(1e10);
  
  obj = optimize(f = ll, interval = c(lower, upper), maximum=TRUE);  # univariate function ll(epsilon) = g(epsilon;y_i)
  epsilon.map = obj$maximum;
  
  ## Compute IL using Laplace approximation
  ## hessian(): Calculate a numerical approximation to the Hessian matrix of a function at a parameter value (here, at epsilon.map)
  ## d2l: g''(epsilon;y_i), calculated numerically --> hessian matrix may not be well-estimated!
  
  d2l = -hessian(ll, epsilon.map);
  il = exp(obj$objective) / sqrt(prior.var * d2l);
  
  list(l = il, epsilon.map = epsilon.map, log.phi.map = prior.log.phi + epsilon.map);
}