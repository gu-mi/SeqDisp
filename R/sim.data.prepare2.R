
#' @title Preparing data for the power-robustness simulations
#' 
#' @details Currently we can simulate three types of dispersion models: "Linear" (NBP), "Quadratic" (NBQ) or "Non-parametric" (Trended).
#' This function differs from \code{sim.data.prepare} by using fold changes from real data estimates instead of specifying a scalar.
#' 
#' @param dt  the full dataset to be loaded and used as a template for the simulation
#' 
#' @param sub.col  column indices used to subset the full dataset
#' 
#' @param grp.ids  group ids to distinguish between the two groups
#' 
#' @param disp  the dispersion type to be simulated. Specify "Quadratic" by default
#' 
#' @param method  estimation method for dispersion, either "ML" or "MAPL" (default is "MAPL")
#' 
#' @param init.fit  (logical) whether to estimate the dispersion model and save as a template (for later use)
#' 
#' @param m  number of genes (default = 5000 genes) subsetted for actual analysis
#' 
#' @param n  number of samples (default = 6 samples)
#' 
#' @param s  library size (default = 1e7, 10 million reads per library)
#' 
#' @param perc.de  percentage of DE genes (default = 0.2)
#' 
#' @param sigma  random normal(0, sigma) noise to phi: phi * exp(rnorm(m, 0, sigma)). 
#' If sigma = 0 (default), data without noise are simulated
#' 
#' @param seed  random generator seed for reproducible research
#' 
#' @param out.path  directory path where the intermediate outputs are stored (and then fetched)
#' 
#' @return a list of read counts ("y"), design matrix ("x"), library size ("s"), indices of DE genes ("idx.DE"),
#' percentage of DE genes ("perc.de"), number of genes ("m"), number of samples ("n"), dispersion ("phi"), and estimated relative 
#' mean frequency ("pi")
#' 
#' @export
#' 
#' @author Gu Mi
#' 
sim.data.prepare2 = function(dt=NULL,
                             sub.col=1:6,
                             grp.ids=c(rep(1,length(sub.col)/2), rep(2,length(sub.col)/2)),
                             disp="Quadratic", 
                             method="MAPL",
                             init.fit=TRUE,
                             m=5000, 
                             n=6, 
                             s=1e7, 
                             perc.de = 0.2, 
                             sigma=0, 
                             seed=539,
                             out.path=getwd()){ 
  
  if (out.path == getwd()) {
    message("You're using current directory to store intermediate outputs. Please specify another output directory if desired.")
  }
  if (is.null(dt)) {
    stop("You must specify an RNA-Seq dataset to be loaded first!")
  }
  if (length(sub.col) != n) {
    stop("The number of columns subsetted does not equal to the number of samples you specified!")
  }
  
  do.call("data", list(dt))
  dt = eval(as.name(dt))
  dt = dt[ ,sub.col]
  dt = as.matrix(dt) # make sure it's a matrix, o.w. Error: INTEGER() can only be applied to a 'integer', not a 'NULL'
  
  # design matrix (two-group comparison case)
  #grp.ids = as.factor(c(rep(1,n/2), rep(2,n/2)))
  grp.ids = as.factor(grp.ids)
  x = model.matrix(~grp.ids)
  
  # filter out very low read counts
  cpms = cpm(dt)   # Returns counts per million (cpm) from a DGEList or matrix object
  keep = rowSums(cpms > 1) >= 1  # exclude very low read counts
  dt2 = dt[keep, ]
  
  # use the entire dataset for the normalization
  nf = estimate.norm.factors(dt2, method="AH2010")
  
  # sub-sample to m genes
  set.seed(seed)
  idx = sample(1:dim(dt2)[1], size=m, replace=FALSE)
  
  # prepare a subset of the entire data
  nb.data = prepare.nb.data(dt2[idx, ], lib.sizes=colSums(dt2), norm.factors=nf)
  
  # for obtaining fold changes from real data, we need to fit a dispersion model first (choose NBQ)
  # estimate NBQ dispersion
  if (init.fit) {
    nbq.dispersion = estimate.dispersion(nb.data, x = x, model="NBQ", method=method)
    save(nbq.dispersion, file=file.path(out.path, "nbq.dispersion.RData"))
  }
  # then, every time, load into R session: nbq.dispersion.RData
  load(file=file.path(out.path, "nbq.dispersion.RData"))
  
  # get beta matrix (m-by-2)
  # note: NBPSeq (version 0.3.0) does not export irls.nb()
  irls.nb.res = NBPSeq:::irls.nb(y=nb.data$counts, s=nb.data$eff.lib.sizes, x=x, phi=nbq.dispersion$estimates, beta0=NA)  
  beta = irls.nb.res$beta
  # THINK: may consider DESeq2's shrinkage approach for LFC estimation
  
  # let non-DE genes have zero beta_2
  n.EE = m * (1-perc.de)   # of EE genes
  set.seed(seed)
  idx.EE = sample(1:m, size=n.EE, replace=FALSE)   # EE gene indices (from the selected m genes)
  beta[idx.EE,2] = 0
  #sum(beta[ ,2] != 0)  # sanity check: = 1000, the number of DE genes if m=5000 and 20% DE
  #sum(beta[ ,2] > 0)   # sanity check: up-regulated DE (about half)
  #sum(beta[ ,2] < 0)   # sanity check: down-regulated DE (about half)
  
  # obtain DE indices used for plotting later
  idx.DE = (1:m)[!((1:m) %in% idx.EE)]
  
  # generate mu
  pi = t(exp(x %*% t(beta)))
  mu = s * pi
  
  # simulate (noised) dispersions according to the QUADRATIC dispersion model pre-estimated
  if (disp == "Quadratic"){
    load(file=file.path(out.path, "nbq.dispersion.RData"))
    z = log(pi[ ,1] / nbq.dispersion$model$offset)   # use first group rel.freq here... a vector, not a matrix
    log.phi.mean = cbind(1, z, z^2) %*% nbq.dispersion$model$par
    phi = exp(log.phi.mean)
    # plot(pi[,1], phi, log="xy")  # sanity check
    if (sigma != 0){
      set.seed(seed)
      phi.noi.vec = phi * exp(rnorm(m, 0, sigma))
      phi = matrix(phi.noi.vec, nr=m, nc=n)
      #plot(pi, phi, log="xy")  # sanity check
      #points(pi[,1], exp(log.phi.mean), col="red")  # true dispersions
    }
  }
  
  # generate responses from Quadratic dispersion model
  set.seed(seed)
  y = rnbinom(m * n, mu=mu, size=1/phi)
  dim(y) = dim(mu)
  rownames(y) = paste0("g", seq(1,m))
  colnames(y) = c( paste0("C", seq(1,sum(grp.ids==1))), paste0("T", seq((sum(grp.ids==1)+1), n)) )
  
  # if a gene has zero counts across all samples, we specify all 1's for that gene (o.w. error from QuasiSeq)
  zero.idx = (rowSums(y) == 0)
  y[zero.idx, ] = 1

  # return quantities as a list
  return(list(counts=y, design=x, lib.sizes = s, idx.DE = idx.DE, perc.de = perc.de, m=m, n=n, phi = phi, pi = pi))
}
