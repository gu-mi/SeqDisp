
##' @title  Simulate datasets to check FDR controls for different NB dispersion models
##' 
##' @description  We simulate read counts of two groups. The (log) mean is determined by normal distribution,
##' and the percentage of DE genes can be specified, with half up- and half down-regulated. Nominal FDR can be
##' specified, and we calculate the (empirical) true FDR and then compare.
##' 
##' @param model  NB dispersion model name
##' @param m  number of genes
##' @param n  number of samples
##' @param mean  the mean of a normal distribution
##' @param sd  the standard deviation of a normal distribution
##' @param perc.de  percentage of DE genes
##' @param FDR  nominal FDR threshold specified below which we call a DE gene
##' @param seed.vec  random number generator seed for reproducibility
##' 
##' @return a data frame with variables to be used for plotting
##' 
##' @examples
##' 
##' library(SeqDisp)
##' library(dplyr)
##' 
##' # check three NB dispersion models with two levels of percentage of DE
##' a1 = sim.fdr.check(model="Genewise", m=5000, n=6, mean=6.5, sd=1.5, perc.de=0.2, FDR=0.1, seed.vec=c(1,2,3))
##' b1 = sim.fdr.check(model="Trended", m=5000, n=6, mean=6.5, sd=1.5, perc.de=0.2, FDR=0.1, seed.vec=c(1,2,3))
##' c1 = sim.fdr.check(model="Tagwise-Trend", m=5000, n=6, mean=6.5, sd=1.5, perc.de=0.2, FDR=0.1, seed.vec=c(1,2,3))
##' a2 = sim.fdr.check(model="Genewise", m=5000, n=6, mean=6.5, sd=1.5, perc.de=0.05, FDR=0.1, seed.vec=c(1,2,3))
##' b2 = sim.fdr.check(model="Trended", m=5000, n=6, mean=6.5, sd=1.5, perc.de=0.05, FDR=0.1, seed.vec=c(1,2,3))
##' c2 = sim.fdr.check(model="Tagwise-Trend", m=5000, n=6, mean=6.5, sd=1.5, perc.de=0.05, FDR=0.1, seed.vec=c(1,2,3))
##' 
##' # combine the results
##' res = rbind(a1,b1,c1,a2,b2,c2)
##' tb_df = tbl_df(data=res)
##' tb_gb = group_by(x=tb_df, model)
##' tb_final = mutate(tb_gb, ratio = actual.fdr/nominal.fdr)
##' 
##' # graphical comparison using ggplot2
##' p = ggplot(data = tb_final, aes(x=model, y=ratio, colour=perc.de, shape=perc.de)) +
##' geom_boxplot(aes(colour=perc.de, shape=perc.de), outlier.size = 0.5, size=0.2) + 
##' geom_hline(yintercept=1, colour="Black", size=0.2, linetype=2) +
##' theme_bw() +
##' scale_x_discrete("") +
##' scale_y_continuous("Actual FDR / Nominal FDR", limits = c(0, max(tb_final[ ,"ratio"]))) +
##' scale_colour_discrete(name =  "% DE genes") +
##' scale_shape_discrete(name = "% DE genes") +
##' theme(plot.title = element_text(face="bold", size=6),
##' # axis labels
##' axis.title.x = element_text(face="bold", size=5),
##' axis.title.y = element_text(face="bold", size=5, angle=90),
##' # axis tick labels
##' axis.text.x = element_text(angle=0, vjust=0.5, size=5), 
##' axis.text.y = element_text(size=5),
##' axis.ticks.length = unit(0.5, "mm"),
##' legend.justification=c(1,1), 
##' legend.position=c(1,1),
##' legend.title = element_text(size=6),
##' legend.key = element_blank(),
##' legend.text = element_text(size=4)
##' )
##' print(p)
##' 
##' @export
##' 
##' @author Gu Mi
##' 
sim.fdr.check = function(model="Genewise", m=5000, n=6, mean=6.5, sd=1.5, perc.de=0.05, FDR=0.1, seed.vec=NULL){
  
  n.sim = length(seed.vec)
  
  # two-group design matrix
  grp.ids = c(rep(1,n/2), rep(2,n/2))
  x = model.matrix(~grp.ids)
  
  # result matrix
  #res.mat = matrix(0, nrow=m*n.sim, ncol=4)
  #colnames(res.mat) = c("p.value","q.value","de.ind","fdr")
  res.mat = matrix(0, nrow=n.sim, ncol=5)
  colnames(res.mat) = c("seed", "model", "perc.de", "nominal.fdr", "actual.fdr")
  res.mat = as.data.frame(res.mat)
  
  for (i in 1:n.sim) {
    
    # mean structure
    set.seed(seed.vec[i])
    mu = exp(rnorm(n=m, mean=mean, sd=sd))
    mu = edgeR:::expandAsMatrix(mu, dim=c(m,n))
    
    # DE genes
    n.de = m * perc.de
    set.seed(seed.vec[i])
    de.idx.up = sample(x=1:m, size=n.de/2, replace=FALSE)
    set.seed(seed.vec[i])
    de.idx.dn = sample(x=(1:m)[!((1:m) %in% de.idx.up)], size=n.de/2, replace=FALSE)
    mu[de.idx.up,((n/2+1):n)] = mu[de.idx.up,((n/2+1):n)] * 2
    mu[de.idx.dn,((n/2+1):n)] = mu[de.idx.dn,((n/2+1):n)] / 2
    de.ind = numeric(m)
    de.ind[c(de.idx.up, de.idx.dn)] = 1
    
    # dispersion phi is determined by 1.5/sqrt(mu)
    phi = 1.5/sqrt(mu)
    
    # simulate read count matrix via mu and phi
    set.seed(seed.vec[i])
    y = rnbinom(n=m*n, size=1/phi, mu=mu)
    dim(y) = c(m,n)
    
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
