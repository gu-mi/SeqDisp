

##' @title  Modified Volcano Plot for RNA-Seq Data Analysis
##' 
##' @description  We use color gradient to show the magnitude of significance.
##' 
##' @param counts  raw RNA-Seq read counts
##' @param x  design matrix
##' @param model  NB dispersion model name (used for testing DE). Specify "NBP" (in the NBPSeq package), "Common", "Genewise", "Trended",
##' "Tagwise-Common" or "Tagwise-Trend" (in the edgeR package)
##' @param title  title of the plot
##' 
##' @return a ggplot2 object
##' 
##' @examples
##' 
##' library(SeqDisp)
##'  
##' data(mouse)
##' 
##' counts = mouse[1:5000,seq(1,6)]  # subset to 5000 genes and two groups
##' grp.ids = as.factor(c(1,1,1,2,2,2))
##' x = model.matrix(~grp.ids)
##' 
##' vocplot.nbp = pl.ma(counts, x, model="NBP")
##' vocplot.gen = pl.ma(counts, x, model="Genewise")
##' 
##' multiplot(vocplot.nbp, vocplot.gen, cols=2, layout=matrix(seq(1,2), nr=1, byrow=TRUE))
##' 
##' @export
##' 
##' @author Gu Mi
##' 

pl.voc = function(counts, x, model=NULL, title=paste0("Volcano Plot: ", model, " Model")){
  
  grp.ids = factor(apply(x, 1, function(x){paste(rev(x), collapse = ".")}), 
                   labels = seq(ncol(x)))
  n.rep.1 = sum(grp.ids == unique(grp.ids)[1])
  first.grp = seq_len(length(grp.ids))[1:n.rep.1]
  second.grp = seq_len(length(grp.ids))[(n.rep.1+1):length(grp.ids)]
  
  expr.dt = cbind(x=apply(counts, 1, function(x) { mean(x[first.grp])}),
                  y=apply(counts, 1, function(x) { mean(x[second.grp])}))
  expr.dt = as.data.frame(expr.dt)
  
  # get significance of each model used
  
  if (model == "NBP") {
    norm.factors = estimate.norm.factors(counts)
    res = nbp.test(counts=counts, grp.ids=grp.ids, grp1=unique(grp.ids)[1], grp2=unique(grp.ids)[2],
                   lib.sizes = colSums(counts), norm.factors = norm.factors, print.level=3)
    expr.dt$logFC = res$log.fc    
    expr.dt$nlog.pv = -log10(res$p.values)
    colnames(expr.dt) = c("first.grp.mean", "second.grp.mean", "logFC", "nlog.pv")  
  }
  if (model == "Common") {
    d = DGEList(counts = counts, group = grp.ids)
    d = calcNormFactors(d, method = "RLE")
    e.com = estimateGLMCommonDisp(d, x)
    com.fit = glmFit(d, x, dispersion = e.com$common.dispersion)
    com.lrt = glmLRT(com.fit, coef = 2)
    expr.dt$logFC = com.lrt$table[ ,"logFC"]
    expr.dt$nlog.pv = -log10(com.lrt$table$PValue)
    colnames(expr.dt) = c("first.grp.mean", "second.grp.mean", "logFC", "nlog.pv")  
  }
  
  if (model == "Genewise") {
    d = DGEList(counts = counts, group = grp.ids)
    d = calcNormFactors(d, method = "RLE")
    e.com = estimateGLMCommonDisp(d, design = x, verbose = FALSE)
    e.gen = estimateGLMTagwiseDisp(d, design = x, dispersion = e.com$common.dispersion, prior.df = 0, trend = FALSE)
    gen.fit = glmFit(d, design = x, dispersion = e.gen$tagwise.dispersion)
    gen.lrt = glmLRT(gen.fit, coef = 2)
    expr.dt$logFC = gen.lrt$table[ ,"logFC"]
    expr.dt$nlog.pv = -log10(gen.lrt$table$PValue)
    colnames(expr.dt) = c("first.grp.mean", "second.grp.mean", "logFC", "nlog.pv")  
  }
  if (model == "Tagwise-Common") {
    d = DGEList(counts = counts, group = grp.ids)
    d = calcNormFactors(d, method = "RLE")
    e.com = estimateGLMCommonDisp(d, x)
    e.tgc = estimateGLMTagwiseDisp(d, design = x, dispersion = e.com$common.dispersion, trend = FALSE)
    tgc.fit = glmFit(d, x, dispersion = e.tgc$tagwise.dispersion)
    tgc.lrt = glmLRT(tgc.fit, coef = 2)
    expr.dt$logFC = tgc.lrt$table[ ,"logFC"]
    expr.dt$nlog.pv = -log10(tgc.lrt$table$PValue)
    colnames(expr.dt) = c("first.grp.mean", "second.grp.mean", "logFC", "nlog.pv")  
  }
  if (model == "Tagwise-Trend") {
    d = DGEList(counts = counts, group = grp.ids)
    d = calcNormFactors(d, method = "RLE")
    e.trd = estimateGLMTrendedDisp(d, design = x)
    e.tgt = estimateGLMTagwiseDisp(d, design = x, dispersion = e.trd$trended.dispersion, trend = TRUE)
    tgt.fit = glmFit(d, design = x, dispersion = e.tgt$tagwise.dispersion)
    tgt.lrt = glmLRT(tgt.fit, coef = 2)
    expr.dt$logFC = tgt.lrt$table[ ,"logFC"]
    expr.dt$nlog.pv = -log10(tgt.lrt$table$PValue)
    colnames(expr.dt) = c("first.grp.mean", "second.grp.mean", "logFC", "nlog.pv")  
  }
  if (model == "Trended") {
    d = DGEList(counts = counts, group = grp.ids)
    d = calcNormFactors(d, method = "RLE")
    e.trd = estimateGLMTrendedDisp(d, design = x)
    trd.fit = glmFit(d, design = x, dispersion = e.trd$trended.dispersion)
    trd.lrt = glmLRT(trd.fit, coef = 2)
    expr.dt$logFC = trd.lrt$table[ ,"logFC"]
    expr.dt$nlog.pv = -log10(trd.lrt$table$PValue)
    colnames(expr.dt) = c("first.grp.mean", "second.grp.mean", "logFC", "nlog.pv")  
  }
  
  # begin volcano plot

  vocplot = 
    ggplot(data=expr.dt, aes(x=logFC, y=nlog.pv)) + 
    # plot genes
    geom_point(data=expr.dt, aes(x=logFC, y=nlog.pv, colour = nlog.pv), alpha=1) +
    # scaling color gradient for significant genes ONLY
    scale_colour_gradient(name="Significance", limits=c(0,5), low = "red", high = "green") +
    scale_y_continuous(limits = c(0, 5)) +
    #
    ggtitle(title) + theme_bw() +
    xlab(expression(log[2]("Fold Change"))) +
    ylab(expression(-log[10]("p-value"))) +
    theme(plot.title = element_text(face = "bold", size = 16), 
          legend.justification = c(0,1), 
          legend.position = c(0,1), 
          legend.title = element_text(face="bold", size=14), 
          legend.key = element_blank(), 
          legend.text = element_text(size = 14))
  
  return(vocplot)
  
}

