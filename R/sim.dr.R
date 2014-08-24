
##' @title Estimate sigma on simulated datasets
##' 
##' @details This function is used for checking the empirical Bayes estimation of dispersion noise sigma on simulated datasets.
##' Two options are available, either to refit the dispersion model being studied, or do not refit (use true dispersions).
##' 
##' @param m  number of genes
##' @param n  number of samples
##' @param x  design matrix
##' @param s  library size
##' @param mu  mean level
##' @param sigma  a vector of true sigma values 
##' @param seed  random number generator seed (for N(0,sigma^2) noise)
##' @param log.phi.mean  a m-by-n matrix of log(phi0) pre-estimated
##' @param model  NB dispersion model name
##' @param method  estimation method for dispersion to be passed to \code{estimate.dispersion.var}, either "ML" or "MAPL" (default is "ML")
##' @param hessian  (logical) whether to calculate the Hessian matrix (default is FALSE)
##' @param refit  (logical) whether to refit the underlying NB dispersion model (default is TRUE)
##' 
##' @return an object from \code{estimate.dispersion.var} (if refit) or \code{optimize.il.prior.var} (if not refit)
##' 
##' @export
##' 
##' @author Yanming Di, Gu Mi
##' 
sim.dr = function(m, n, x, s, mu, sigma, seed, log.phi.mean = NULL, model = NULL, method = "ML", hessian=FALSE, refit = TRUE) {
  
  if (is.null(model)){
    stop("You must specify one of the following models: 'NBP', 'NBQ', 'NBS', 'NB2', 'Trended'!")
  }
  
  if ( refit & is.null(model) ){
    stop("Since you want to refit a model, please specify the model name!")
  }
  
  set.seed(seed)
  ## Add noise (Normal mean 0 and variance sigma^2) --> log.phi.mean = prior.log.phi, is the log(phi0) pre-estimated
  log.phi =  log.phi.mean + rnorm(m, mean = 0, sd = sigma)
  phi = exp(log.phi)
  
  y = rnbinom(m*n, size = 1/phi, mu=mu)
  dim(y) = c(m,n)
  
  ## if we refit the dispersion model ##
  if (refit){
    #
    nb.data = prepare.nb.data(counts=y, lib.sizes=s)
    #
    if (model == "NBP"){
      est.disp = estimate.dispersion(nb.data, x, model="NBP", method=method)
      obj = estimate.dispersion.var(nb.data, dispersion=est.disp, x, hessian=hessian)  
    }
    #
    if (model == "NBQ"){
      est.disp = estimate.dispersion(nb.data, x, model="NBQ", method=method)
      obj = estimate.dispersion.var(nb.data, dispersion=est.disp, x, hessian=hessian)  
    }
    #
    if (model == "NBS"){
      est.disp = estimate.dispersion(nb.data, x, model="NBS", method=method)
      obj = estimate.dispersion.var(nb.data, dispersion=est.disp, x, hessian=hessian)  
    }
    #
    if (model == "NB2"){
      est.disp = estimate.dispersion(nb.data, x, model="NB2", method=method)
      obj = estimate.dispersion.var(nb.data, dispersion=est.disp, x, hessian=hessian)  
    }
    #
    if (model == "Trended"){
      y.dge = DGEList(counts = y, lib.size=s)
      y.dge = estimateGLMTrendedDisp(y.dge, design = x)
      dispersion = expandAsMatrix(x=y.dge$trended.dispersion, dim=c(m,n))
      obj = estimate.dispersion.var.edgeR(nb.data, dispersion=dispersion, x, hessian=hessian)
    }
  }
  
  ## if we do NOT refit the dispersion model ##
  if (!refit) {
    #args(optimize.il.prior.var)
    obj = optimize.il.prior.var(y, s, x, log.phi.mean, hessian=hessian)  # pass log.phi.mean as the log.phi argument it's "prior.log.phi"
    obj  # $maximum: gives the sigma.hat^2
  }
  
  return(obj)
}


