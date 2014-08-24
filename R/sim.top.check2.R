

##' @title  Simulate datasets (from real data) to check top DE genes for different NB dispersion models
##' 
##' @description  We simulate read counts of two groups from a real dataset. The (log) mean is determined by normal distribution,
##' and the percentage of DE genes can be specified, with half up- and half down-regulated. The number of DE genes among the 
##' pre-specified top genes (default 500) is reported.
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
##' @param top  the number of top genes (from DE test p-values) considered
##' 
##' @return a data frame with variables to be used for plotting
##' 
##' @export
##' 
##' @author Gu Mi
##' 
sim.top.check2 = function (model = "Genewise", dt = "mouse", sub.col = 1:6, disp.vec = "Quadratic", 
                           init.fit = FALSE, m = 5000, n = 6, s = 1.5e+07, perc.de = NULL, 
                           fc.vec = NULL, noise.vec = 1.05, seed.vec = NULL, out.path = path.o, 
                           top = 500) 
{
  n.sim = length(seed.vec)
  res.mat = matrix(0, nrow = n.sim, ncol = 4)
  colnames(res.mat) = c("seed", "model", "perc.de", "n.true.de")
  res.mat = as.data.frame(res.mat)
  for (i in 1:n.sim) {
    d.sim = sim.data.prepare2(dt = dt, sub.col = sub.col, 
                              disp = disp.vec, init.fit = init.fit, m = m, n = n, 
                              s = s, perc.de = perc.de, sigma = noise.vec, seed = seed.vec[i], 
                              out.path = path.o)
    y = d.sim$counts
    # to avoid fitting issues in QuasiSeq, we add 1 to all counts
    # error message in QuasiSeq: 1 genes have 0 counts across all samples. Please remove genes with zero total counts before analyzing.
    y = y + 1
    x = d.sim$design
    DE.vec = numeric(m)
    DE.vec[d.sim$idx.DE] = 1
    
    grp.ids = factor(apply(x, 1, function(x) {
      paste(rev(x), collapse = ".")
    }), labels = seq(ncol(x)))
    de.ind = numeric(m)
    de.ind[d.sim$idx.DE] = 1
    if (model == "NBQ") {
      nbq.fit = nb.glm.test(counts = y, x = x, beta0 = c(NA, 0), normalization.method = "AH2010", dispersion.model = "NBQ", 
                            tests = "HOA")
      pval.y = nbq.fit$test.results$HOA$p.values
      o.pv = order(pval.y)
      n.true.de = sum(DE.vec[o.pv][1:top], na.rm=TRUE)
    }
    if (model == "Common") {
      d = DGEList(counts = y, group = grp.ids)
      d = calcNormFactors(d, method = "RLE")
      e.com = estimateGLMCommonDisp(d, x)
      com.fit = glmFit(d, x, dispersion = e.com$common.dispersion)
      com.lrt = glmLRT(com.fit, coef = 2)
      pval.y = com.lrt$table$PValue
      o.pv = order(pval.y)
      n.true.de = sum(DE.vec[o.pv][1:top], na.rm=TRUE)
    }
    if (model == "Genewise") {
      d = DGEList(counts = y, group = grp.ids)
      d = calcNormFactors(d, method = "RLE")
      e.com = estimateGLMCommonDisp(d, design = x, verbose = FALSE)
      e.gen = estimateGLMTagwiseDisp(d, design = x, dispersion = e.com$common.dispersion, 
                                     prior.df = 0, trend = FALSE)
      gen.fit = glmFit(d, design = x, dispersion = e.gen$tagwise.dispersion)
      gen.lrt = glmLRT(gen.fit, coef = 2)
      pval.y = gen.lrt$table$PValue
      o.pv = order(pval.y)
      n.true.de = sum(DE.vec[o.pv][1:top], na.rm=TRUE)
    }
    if (model == "Trended") {
      d = DGEList(counts = y, group = grp.ids)
      d = calcNormFactors(d, method = "RLE")
      e.trd = estimateGLMTrendedDisp(d, design = x)
      trd.fit = glmFit(d, design = x, dispersion = e.trd$trended.dispersion)
      trd.lrt = glmLRT(trd.fit, coef = 2)
      pval.y = trd.lrt$table$PValue
      o.pv = order(pval.y)
      n.true.de = sum(DE.vec[o.pv][1:top], na.rm=TRUE)
    }
    if (model == "Tagwise-Trend") {
      d = DGEList(counts = y, group = grp.ids)
      d = calcNormFactors(d, method = "RLE")
      d = estimateGLMTrendedDisp(d, design = x)
      e.tgt = estimateGLMTagwiseDisp(d, design = x, trend = TRUE)
      tgt.fit = glmFit(d, design = x, dispersion = e.tgt$tagwise.dispersion)
      tgt.lrt = glmLRT(tgt.fit, coef = 2)
      pval.y = tgt.lrt$table$PValue
      o.pv = order(pval.y)
      n.true.de = sum(DE.vec[o.pv][1:top], na.rm=TRUE)
    }
    if (model == "QL") {
      design.list = vector("list", 2)
      design.list[[1]] = x
      design.list[[2]] = rep(1, length(grp.ids))
      log.offset = log(apply(y, 2, quantile, 0.75))
      fit = QL.fit(counts = y, design.list = design.list, 
                   log.offset = log.offset, Model = "NegBin", method = "glm", 
                   NBdisp = "trend", print.progress = FALSE)
      res = QL.results(fit, Plot = FALSE)
      pval.y = as.numeric(res$P.values$QL)
      o.pv = order(pval.y)
      n.true.de = sum(DE.vec[o.pv][1:top], na.rm=TRUE)
    }
    if (model == "QLShrink") {
      design.list = vector("list", 2)
      design.list[[1]] = x
      design.list[[2]] = rep(1, length(grp.ids))
      log.offset = log(apply(y, 2, quantile, 0.75))
      fit = QL.fit(counts = y, design.list = design.list, 
                   log.offset = log.offset, Model = "NegBin", method = "glm", 
                   NBdisp = "trend", print.progress = FALSE)
      res = QL.results(fit, Plot = FALSE)
      pval.y = as.numeric(res$P.values$QLShrink)
      o.pv = order(pval.y)
      n.true.de = sum(DE.vec[o.pv][1:top], na.rm=TRUE)
    }
    if (model == "QLSpline") {
      design.list = vector("list", 2)
      design.list[[1]] = x
      design.list[[2]] = rep(1, length(grp.ids))
      log.offset = log(apply(y, 2, quantile, 0.75))
      fit = QL.fit(counts = y, design.list = design.list, 
                   log.offset = log.offset, Model = "NegBin", method = "glm", 
                   NBdisp = "trend", print.progress = FALSE)
      res = QL.results(fit, Plot = FALSE)
      pval.y = as.numeric(res$P.values$QLSpline)
      o.pv = order(pval.y)
      n.true.de = sum(DE.vec[o.pv][1:top], na.rm=TRUE)
    }

    res.mat[i, "seed"] = seed.vec[i]
    res.mat[i, "model"] = model
    res.mat[i, "perc.de"] = as.character(perc.de)
    res.mat[i, "n.true.de"] = n.true.de
  }
  return(res.mat)
}
