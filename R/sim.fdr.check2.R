##' @title  Simulate datasets (from real data) to check FDR controls for different NB dispersion models
##' 
##' @description  We simulate read counts of two groups from a real dataset. The (log) mean is determined by normal distribution,
##' and the percentage of DE genes can be specified, with half up- and half down-regulated. Nominal FDR can be
##' specified, and we calculate the (empirical) true FDR and then compare.
##' 
##' @param model  NB dispersion model name
##' @param dt  data name
##' @param sub.col  subset of columns
##' @param disp.vec  dispersion name
##' @param init.fit  whether to have an initial fit for the dispersion model
##' @param m  number of genes
##' @param n  number of samples
##' @param s  library size
##' @param perc.de  percentage of DE genes
##' @param fc.vec  fold change
##' @param noise.vec  noise
##' @param seed.vec  random number generator seed for reproducibility
##' @param out.path  output directory
##' @param FDR  nominal FDR threshold specified below which we call a DE gene
##' 
##' @return a data frame with variables to be used for plotting
##' 
##' @export
##' 
##' @author Gu Mi
##' 
sim.fdr.check2 = function(model="Genewise", 
                          dt = "mouse",
                          sub.col = 1:6,
                          disp.vec = "Quadratic",
                          init.fit = FALSE,
                          m = 5000,
                          n = 6,
                          s = 1.5e7, 
                          perc.de = NULL,
                          fc.vec = NULL,
                          noise.vec = 1.05,
                          seed.vec=NULL,
                          out.path = path.o,
                          FDR=0.1) {
  
  n.sim = length(seed.vec)
  
  # result matrix
  res.mat = matrix(0, nrow=n.sim, ncol=5)
  colnames(res.mat) = c("seed", "model", "perc.de", "nominal.fdr", "actual.fdr")
  res.mat = as.data.frame(res.mat)
  
  for (i in 1:n.sim) {
    
    d.sim = sim.data.prepare2(dt = dt, sub.col = sub.col, disp = disp.vec, init.fit = init.fit, m = m, n = n, 
                              s = s, perc.de = perc.de, sigma = noise.vec, seed = seed.vec[i], out.path = path.o) 
    
    # extract information and quantities
    y = d.sim$counts
    x = d.sim$design
    grp.ids = factor(apply(x, 1, function(x){paste(rev(x), collapse = ".")}), labels = seq(ncol(x)))
    de.ind = numeric(m)
    de.ind[d.sim$idx.DE] = 1
    
    if (model == "NBQ") {
      nbq.fit = nb.glm.test(counts = y, x = x, beta0 = c(NA, 0), normalization.method = "AH2010", dispersion.model = "NBQ", tests = "HOA")
      pval.y = nbq.fit$test.results$HOA$p.values
      qval.y = p.adjust(p=pval.y, method="BH")
    }
    
    if (model == "Common") {
      d = DGEList(counts = y, group = grp.ids)
      d = calcNormFactors(d, method = "RLE")
      e.com = estimateGLMCommonDisp(d, x)
      com.fit = glmFit(d, x, dispersion = e.com$common.dispersion)
      com.lrt = glmLRT(com.fit, coef = 2)
      pval.y = com.lrt$table$PValue
      qval.y = p.adjust(p=pval.y, method="BH")
    }
    
    if (model == "Genewise"){
      # genewise tests for DE genes
      d = DGEList(counts = y, group = grp.ids)
      d = calcNormFactors(d, method="RLE")
      e.com = estimateGLMCommonDisp(d, design=x, verbose=FALSE) 
      e.gen = estimateGLMTagwiseDisp(d, design=x, dispersion=e.com$common.dispersion, prior.df=0, trend=FALSE)  # prior.df = 0
      # e.gen: trend is FALSE since we didn't use estimateGLMTrendedDisp() beforehand
      gen.fit = glmFit(d, design=x, dispersion=e.gen$tagwise.dispersion)
      gen.lrt = glmLRT(gen.fit,coef=2)
      pval.y = gen.lrt$table$PValue
      qval.y = p.adjust(p=pval.y, method="BH")
    }
    
    if (model == "Trended"){
      d = DGEList(counts = y, group = grp.ids)
      d = calcNormFactors(d, method = "RLE")
      e.trd = estimateGLMTrendedDisp(d, design = x)
      trd.fit = glmFit(d, design = x, dispersion = e.trd$trended.dispersion)
      trd.lrt = glmLRT(trd.fit, coef = 2)
      pval.y = trd.lrt$table$PValue
      qval.y = p.adjust(p=pval.y, method="BH")
    }
    
    if (model == "Tagwise-Trend") {
      d = DGEList(counts = y, group = grp.ids)
      d = calcNormFactors(d, method = "RLE")
      e.trd = estimateGLMTrendedDisp(d, design = x)
      e.tgt = estimateGLMTagwiseDisp(d, design = x, dispersion = e.trd$trended.dispersion, trend = TRUE)
      tgt.fit = glmFit(d, design = x, dispersion = e.tgt$tagwise.dispersion)
      tgt.lrt = glmLRT(tgt.fit, coef = 2)
      pval.y = tgt.lrt$table$PValue
      qval.y = p.adjust(p=pval.y, method="BH")
    }
    
    if (model == "QL"){
      design.list = vector("list",2)
      design.list[[1]] = x 
      design.list[[2]] = rep(1,length(grp.ids))
      log.offset = log(apply(y,2,quantile,0.75))
      fit = QL.fit(counts=y, design.list=design.list, log.offset=log.offset, Model="NegBin",method="glm", NBdisp="trend", print.progress=FALSE)
      res = QL.results(fit, Plot=FALSE)
      pval.y = as.numeric(res$P.values$QL)
      qval.y = p.adjust(p=pval.y, method="BH")
    }
    
    if (model == "QLShrink"){
      design.list = vector("list",2)
      design.list[[1]] = x  
      design.list[[2]] = rep(1,length(grp.ids))
      log.offset = log(apply(y,2,quantile,0.75))
      fit = QL.fit(counts=y, design.list=design.list, log.offset=log.offset, Model="NegBin",method="glm", NBdisp="trend", print.progress=FALSE)
      res = QL.results(fit, Plot=FALSE)
      pval.y = as.numeric(res$P.values$QLShrink)
      qval.y = p.adjust(p=pval.y, method="BH")
    }
    
    if (model == "QLSpline"){
      design.list = vector("list",2)
      design.list[[1]] = x  
      design.list[[2]] = rep(1,length(grp.ids))
      log.offset = log(apply(y,2,quantile,0.75))
      fit = QL.fit(counts=y, design.list=design.list, log.offset=log.offset, Model="NegBin",method="glm", NBdisp="trend", print.progress=FALSE)
      res = QL.results(fit, Plot=FALSE)
      pval.y = as.numeric(res$P.values$QLSpline)
      qval.y = p.adjust(p=pval.y, method="BH")
    }
    
    fp = sum(qval.y < FDR & de.ind == 0)
    tp = sum(qval.y < FDR & de.ind == 1)
    fdr = fp / (tp + fp) 
    
    # fill out the result matrix
    res.mat[i,"seed"] = seed.vec[i]
    res.mat[i,"model"] = model
    res.mat[i,"perc.de"] = as.character(perc.de)
    res.mat[i,"nominal.fdr"] = FDR
    res.mat[i,"actual.fdr"] = fdr
  }
  
  return(res.mat)
}
