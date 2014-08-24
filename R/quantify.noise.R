
##' @title Prepare NB dispersions for NBQ dispersion models, and calculate "sigma" in real and simulated datasets
##' 
##' @details The dispersion results are saved as RData files, which will be used in subsequent steps to calculate
##' 
##'   log.phi.mean = cbind(1, z, z^2) %*% dispersion$model$par
##'   
##' e.g. the "dispersion" object
##' 
##' @param dt  name of the dataset
##' @param sigma.v  a vector of true sigma
##' @param m  number of genes to subset
##' @param sub.col  indices for subsetting the columns of dt
##' @param grp.ids  group ids to distinguish between the two groups
##' @param evaluate  a vector of logicals for evaluating underlying dispersion, real data sigma, and simulated data sigma's
##' @param method  estimation method for dispersion to be used for quantify residual dispersion variation, either "ML" or "MAPL"
##' (default is "ML")
##' @param hessian  whether to calculate the hessian matrix for getting standard error of the estimated sigma (default is FALSE)
##' @param seed  random generator seed for sampling a subset of genes
##' @param path.o  output directory
##' 
##' @return many saved intermediate files used for quantifying noise
##' 
##' @author Gu Mi
##' 
##' @export
##' 
quantify.noise = function(dt=NULL, 
                          sigma.v = c(0.5, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.5),
                          m=5000, 
                          sub.col=NULL, 
                          grp.ids=c(rep(1,length(sub.col)/2), rep(2,length(sub.col)/2)),
                          evaluate=c(TRUE, FALSE, FALSE),
                          method="ML",
                          hessian=FALSE,
                          seed=539, 
                          path.o=getwd()) {
  
  setwd(path.o)
  
  do.call("data", list(dt))
  dta = eval(as.name(dt))
  dta = dta[, sub.col]
  
  n = dim(dta)[2]
  #grp.ids = as.factor(c(rep(1,n/2), rep(2,n/2)))
  grp.ids = as.factor(grp.ids)
  x = model.matrix(~grp.ids)
  
  # we estimate the dispersion noise sigma for the control AND treatment groups to a subset of genes
  cpms = cpm(dta)
  keep = rowSums(cpms > 1) >= 1
  tmp = dta[keep, ]
  set.seed(seed)
  idx = sample(1:dim(tmp)[1], size=m, replace=FALSE)
  dta.sub = as.matrix(tmp[idx, ])
  
  # use the FULL dataset to calculate normalization factor "nf"
  nf = estimate.norm.factors(tmp, method="AH2010")
  # this should be the same (in edgeR):
  # calcNormFactors(tmp, method="RLE")
  
  # prepare nb data: use subset of 5000 genes in the first argument, but pass "nf" to norm.factors
  nb.data = prepare.nb.data(dta.sub, lib.sizes=colSums(tmp), norm.factors=nf)
  
  if (evaluate[1]){
    # assume NBQ as the underlying dispersion model
    dispersion = estimate.dispersion(nb.data=nb.data, x=x, model="NBQ", method=method)
    save(dispersion, file=file.path(path.o, paste0(dt, ".", m, ".NBQ.RData")))   # e.g. mouse.5000.NBQ.RData
  }
  
  load(file=file.path(path.o, paste0(dt, ".", m, ".NBQ.RData")))
  
  if (evaluate[2]){
    # estimate sigma for real data
    obj = estimate.dispersion.var(nb.data=nb.data, dispersion=dispersion, x=x, hessian=hessian)
    # save object obj
    file.out = sprintf(file.path(path.o, "sigma.m%d.n%d.NBQ.real.Rdata"), m, n)
    save(obj, file=file.out)
    
    # save estimated sigma from real data
    sigma.real = sqrt(obj$maximum) 
    real.data.sigma.obj = data.frame(sigma.v = sigma.v, real.data.sigmas = rep(sigma.real, length(sigma.v)))
    write.table(real.data.sigma.obj, file=file.path(path.o, paste0(dt, ".sigma.obj.NBQ.real.txt")), 
                quote=FALSE, row.names=FALSE, col.names=TRUE)
    if (hessian) {
      se = sqrt(- 1/obj$hessian)
      write.table(se, file=file.path(path.o, paste0(dt, "se.sigma.NBQ.real.txt")), quote=FALSE, row.names=FALSE, col.names=FALSE)
    }
  }
  
  if (evaluate[3]) {
    #
    # Simulated datasets using one of the models (repeat 3 times: plot the median value in blue and others in gray)
    #
    # basic info. on real data
    mu.dta = mean(log(rowMeans(tmp)))
    sd.dta = sd(log(rowMeans(tmp)))
    s.dta = max(colSums(tmp))
    
    set.seed(seed)
    mu = sort(rlnorm(m, mu.dta, sd.dta))
    s = rep(s.dta, n)
    pi = mu / s[1]
    pi = edgeR:::expandAsMatrix(x=pi, dim=c(m,n))
    
    # Simulated datasets #
    # prepare for the estimated sigma's matrix (with meta-data)
    n.sigma = length(sigma.v)
    n.rep = 3  # at each true sigma, we evaluate 3 times
    n.sim = n.sigma * n.rep
    res = matrix(NA, n.sim, 3)
    colnames(res) = c("seed", "sigma", "sigma.hat")
    res[ ,"seed"] = 5000 + (1:n.sim) * 5000
    res[ ,"sigma"] = rep(sigma.v, each=n.rep)
    
    # prepare for log.phi.mean for different dispersion models (now only include NBQ)
    #z = log(pi[ ,1] / dispersion$model$pi.offset);
    # NBPSeq v.0.3.5: pi.offset --> offset
    z = log(pi[ ,1] / dispersion$model$offset);
    log.phi.mean = cbind(1, z, z^2) %*% dispersion$model$par  # log.phi.mean = prior.log.phi, is the log(phi0) pre-estimated
    #plot(pi[,1], exp(log.phi.mean), log="xy")  # sanity check
    
    # begin sim.dr for estimating sigma for a chosen model
    for (i in 1:n.sim) {
      obj = sim.dr(m=m, n=n, x=x, s=s, mu=mu, sigma=res[i,"sigma"], seed=res[i,"seed"], log.phi.mean=log.phi.mean, model="NBQ", 
                   method=method, hessian=hessian, refit=TRUE)
      res[i,3] = sqrt(obj$maximum)
      print(res[i, ])
    }  
    file.out = sprintf(file.path(path.o, "sigma.m%d.n%d.NBQ.sim.RData"), m, n)
    save(res, file=file.out)
    
    sim.sigma.res = as.data.frame(res)  
    
    sim.sigma.v.medians = numeric(length(sigma.v))
    for (i in seq_len(length(sigma.v))){
      res_df = tbl_df(data=sim.sigma.res)
      res_gb = group_by(x=res_df, sigma)
      sim.sigma.v.medians[i] = as.data.frame(summarise(res_gb, median_sigma = median(sigma.hat, na.rm=TRUE)))$median_sigma[i]
    }
    sim.sigma.obj = data.frame(rep(sigma.v, n.rep), rep(sim.sigma.v.medians, n.rep), sim.sigma.res)
    colnames(sim.sigma.obj) = c("sigma.v", "sim.sigma.v.medians", "seed", "sigma", "sigma.hat")
    
    write.table(sim.sigma.obj, file=file.path(path.o, paste0(dt, ".sigma.obj.NBQ.sim.txt")), 
                quote=FALSE, row.names=FALSE, col.names=TRUE)
  }
  
}
