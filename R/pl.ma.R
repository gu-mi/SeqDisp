
##' @title  Modified MA Plot for RNA-Seq Data Analysis
##' 
##' @description  Typical MA plot would not have the significance. This modified function introduces the significance 
##' used in Volcano plot into MA plot for the selected DE genes (by FDR), to better reveal the relationship between test significance 
##' and gene expression levels (fold change). We use color gradient to show the magnitude of FDR (red: most significant;
##' blue: significant but close to the FDR threshold). In addition, for the genewise approach, we use color gradient to show the dispersion estimates for the top selected DE genes. For other approaches, we plot in red the top-selected genes (by the \code{top} argument).
##' 
##' @param counts  raw RNA-Seq read counts
##' @param x  design matrix
##' @param model  NB dispersion model name (used for testing DE). Specify "NBP", "NBQ", "NB2", "Genewise-HOA" (in the NBPSeq package),
##' "Common", "Genewise", "Trended", "Tagwise-Common" or "Tagwise-Trend" (in the edgeR package)
##' @param FDR  threshold of false discovery rate (default = 0.1) to call a gene significant
##' @param top  a scalar showing the number of most significant DE genes
##' @param gradient  (logical) whether to use color gradient to show FDR for DE genes (if \code{FDR=TRUE}), or show estimated
##' dispersions for the genewise model (if \code{top=TRUE}). If \code{FALSE}, DE genes will be shown in solid red.
##' @param title  title of the plot
##' 
##' @return a ggplot2 object
##' 
##' @examples
##' 
##' library(SeqDisp)
##'  
##' data(mouse)
##' counts = mouse[1:5000,seq(1,6)]  # subset to 5000 genes and two groups
##' grp.ids = as.factor(c(1,1,1,2,2,2))
##' x = model.matrix(~grp.ids)
##' 
##' # select DE genes by FDR threshold
##' map.nbp.fdr = pl.ma(counts, x, model="NBP", FDR=0.1, gradient=TRUE)
##' map.gen.fdr = pl.ma(counts, x, model="Genewise", FDR=0.1, gradient=TRUE)
##' multiplot(map.nbp.fdr, map.gen.fdr, cols=2, layout=matrix(seq(1,2), nr=1, byrow=TRUE))
##' 
##' # select top DE genes
##' map.nbp.top = pl.ma(counts, x, model="NBP", top=100)
##' map.gen.top = pl.ma(counts, x, model="Genewise", top=100)
##' multiplot(map.nbp.top, map.gen.top, cols=2, layout=matrix(seq(1,2), nr=1, byrow=TRUE))
##' 
##' @export
##' 
##' @author Gu Mi
##' 

pl.ma = function(counts, x, model=NULL, FDR=NULL, top=NULL, gradient=FALSE, title=paste0("MA Plot: ", model, " Model")){
  
  m = dim(counts)[1]
  n = dim(counts)[2]
  
  if (all(is.null(c(top, FDR)))) {
    stop("Need value for top genes (via argument 'top') or FDR cut-off (via argument 'FDR'.")
  }
  if (!is.null(top) & !is.null(FDR)) {
    stop("Please specify a value to either 'FDR' or 'top', not both.")
  }
  
  grp.ids = factor(apply(x, 1, function(x){paste(rev(x), collapse = ".")}), 
                   labels = seq(ncol(x)))
  n.rep.1 = sum(grp.ids == unique(grp.ids)[1])
  first.grp = seq_len(length(grp.ids))[1:n.rep.1]
  second.grp = seq_len(length(grp.ids))[(n.rep.1+1):length(grp.ids)]
  
  expr.dt = cbind(first.grp.mean=apply(counts, 1, function(x) { mean(x[first.grp])}),
                  second.grp.mean=apply(counts, 1, function(x) { mean(x[second.grp])}))
  expr.dt = as.data.frame(expr.dt)
  
  # get significance of each model used and dispersion estimates (for genewise approach)
  
  if (model %in% c("NBP","NBQ","NB2")) {
    nf = estimate.norm.factors(counts=counts, lib.sizes=colSums(counts), method="AH2010")
    prep.nbp = prepare.nbp(counts=counts, grp.ids=grp.ids, lib.sizes=colSums(counts), norm.factors=nf, print.level=3)
    e.disp = estimate.disp(prep.nbp, model = model, print.level = 3)
    res = exact.nb.test(e.disp, grp1=unique(grp.ids)[1], grp2=unique(grp.ids)[2], print.level = 3)
    expr.dt$logFC = res$log.fc
    expr.dt$logCPM = 1e+06 * res$pooled.pie
    expr.dt$FDR = res$q.values
    expr.dt = expr.dt[complete.cases(expr.dt), ]
    # if ask for DE genes by FDR (to be highlighted)
    if (!is.null(FDR)){
      expr.dt$sig.ind = ifelse(test=expr.dt$FDR < FDR, 1, 0)
      non.de.dt = expr.dt[expr.dt$sig.ind == 0, ]
      de.dt = expr.dt[expr.dt$sig.ind == 1, ]
    }
    # if ask for top DE genes (to be highlighted)
    if (!is.null(top)){
      expr.dt$top.ind = numeric(dim(expr.dt)[1])
      pv.idx = order(expr.dt$FDR)
      expr.dt = expr.dt[pv.idx, ]
      expr.dt[(1:top),"top.ind"] = 1
      non.top.dt = expr.dt[expr.dt$top.ind == 0, ]
      top.dt = expr.dt[expr.dt$top.ind == 1, ]
    }
  }
  
  if (model == "Common") {
    d = DGEList(counts = counts, group = grp.ids)
    d = calcNormFactors(d, method = "RLE")
    e.com = estimateGLMCommonDisp(d, x)
    com.fit = glmFit(d, x, dispersion = e.com$common.dispersion)
    com.lrt = glmLRT(com.fit, coef = 2)
    expr.dt$logFC = com.lrt$table[ ,"logFC"]
    expr.dt$logCPM = com.lrt$table[ ,"logCPM"]
    expr.dt$FDR = p.adjust(p=com.lrt$table$PValue, method="BH")
    expr.dt = expr.dt[complete.cases(expr.dt), ]
    # if ask for DE genes by FDR (to be highlighted)
    if (!is.null(FDR)){
      expr.dt$sig.ind = ifelse(test=expr.dt$FDR < FDR, 1, 0)
      non.de.dt = expr.dt[expr.dt$sig.ind == 0, ]
      de.dt = expr.dt[expr.dt$sig.ind == 1, ]
    }
    # if ask for top DE genes (to be highlighted)
    if (!is.null(top)){
      expr.dt$top.ind = numeric(dim(expr.dt)[1])
      pv.idx = order(com.lrt$table$PValue)
      expr.dt = expr.dt[pv.idx, ]
      expr.dt[(1:top),"top.ind"] = 1
      non.top.dt = expr.dt[expr.dt$top.ind == 0, ]
      top.dt = expr.dt[expr.dt$top.ind == 1, ]
    }
  }
  
  if (model == "Genewise") {
    d = DGEList(counts = counts, group = grp.ids)
    d = calcNormFactors(d, method = "RLE")
    e.com = estimateGLMCommonDisp(d, design = x, verbose = FALSE)
    e.gen = estimateGLMTagwiseDisp(d, design = x, dispersion = e.com$common.dispersion, prior.df = 0, trend = FALSE)
    gen.fit = glmFit(d, design = x, dispersion = e.gen$tagwise.dispersion)
    gen.lrt = glmLRT(gen.fit, coef = 2)
    expr.dt$logFC = gen.lrt$table[ ,"logFC"]
    expr.dt$logCPM = gen.lrt$table[ ,"logCPM"]
    expr.dt$FDR = p.adjust(p=gen.lrt$table$PValue, method="BH")
    expr.dt = expr.dt[complete.cases(expr.dt), ]
    # if ask for DE genes by FDR (to be highlighted)
    if (!is.null(FDR)){
      expr.dt$sig.ind = ifelse(test=expr.dt$FDR < FDR, 1, 0)
      non.de.dt = expr.dt[expr.dt$sig.ind == 0, ]
      de.dt = expr.dt[expr.dt$sig.ind == 1, ]
    }
    # if ask for top DE genes (to be highlighted)
    if (!is.null(top)){
      expr.dt$top.ind = numeric(dim(expr.dt)[1])
      pv.idx = order(gen.lrt$table$PValue)
      expr.dt = expr.dt[pv.idx, ]
      expr.dt[(1:top),"top.ind"] = 1
      expr.dt$est.disp = e.gen$tagwise.dispersion
      non.top.dt = expr.dt[expr.dt$top.ind == 0, ]
      top.dt = expr.dt[expr.dt$top.ind == 1, ]
    }
  }
  
  if (model == "Tagwise-Common") {
    d = DGEList(counts = counts, group = grp.ids)
    d = calcNormFactors(d, method = "RLE")
    e.com = estimateGLMCommonDisp(d, x)
    e.tgc = estimateGLMTagwiseDisp(d, design = x, dispersion = e.com$common.dispersion, trend = FALSE)
    tgc.fit = glmFit(d, x, dispersion = e.tgc$tagwise.dispersion)
    tgc.lrt = glmLRT(tgc.fit, coef = 2)
    expr.dt$logFC = tgc.lrt$table[ ,"logFC"]
    expr.dt$logCPM = tgc.lrt$table[ ,"logCPM"]
    expr.dt$FDR = p.adjust(p=tgc.lrt$table$PValue, method="BH")
    expr.dt = expr.dt[complete.cases(expr.dt), ]
    # if ask for DE genes by FDR (to be highlighted)
    if (!is.null(FDR)){
      expr.dt$sig.ind = ifelse(test=expr.dt$FDR < FDR, 1, 0)
      colnames(expr.dt) = c("first.grp.mean", "second.grp.mean", "logFC", "logCPM", "FDR", "sig.ind")  
      non.de.dt = expr.dt[expr.dt$sig.ind == 0, ]
      de.dt = expr.dt[expr.dt$sig.ind == 1, ]
    }
    # if ask for top DE genes (to be highlighted)
    if (!is.null(top)){
      expr.dt$top.ind = numeric(dim(expr.dt)[1])
      pv.idx = order(tgc.lrt$table$PValue)
      expr.dt = expr.dt[pv.idx, ]
      expr.dt[(1:top),"top.ind"] = 1
      non.top.dt = expr.dt[expr.dt$top.ind == 0, ]
      top.dt = expr.dt[expr.dt$top.ind == 1, ]
    }
  }
  
  if (model == "Tagwise-Trend") {
    d = DGEList(counts = counts, group = grp.ids)
    d = calcNormFactors(d, method = "RLE")
    e.trd = estimateGLMTrendedDisp(d, design = x)
    e.tgt = estimateGLMTagwiseDisp(d, design = x, dispersion = e.trd$trended.dispersion, trend = TRUE)
    tgt.fit = glmFit(d, design = x, dispersion = e.tgt$tagwise.dispersion)
    tgt.lrt = glmLRT(tgt.fit, coef = 2)
    expr.dt$logFC = tgt.lrt$table[ ,"logFC"]
    expr.dt$logCPM = tgt.lrt$table[ ,"logCPM"]
    expr.dt$FDR = p.adjust(p=tgt.lrt$table$PValue, method="BH")
    expr.dt = expr.dt[complete.cases(expr.dt), ]
    # if ask for DE genes by FDR (to be highlighted)
    if (!is.null(FDR)){
      expr.dt$sig.ind = ifelse(test=expr.dt$FDR < FDR, 1, 0)
      colnames(expr.dt) = c("first.grp.mean", "second.grp.mean", "logFC", "logCPM", "FDR", "sig.ind")  
      non.de.dt = expr.dt[expr.dt$sig.ind == 0, ]
      de.dt = expr.dt[expr.dt$sig.ind == 1, ]
    }
    # if ask for top DE genes (to be highlighted)
    if (!is.null(top)){
      expr.dt$top.ind = numeric(dim(expr.dt)[1])
      pv.idx = order(tgt.lrt$table$PValue)
      expr.dt = expr.dt[pv.idx, ]
      expr.dt[(1:top),"top.ind"] = 1
      non.top.dt = expr.dt[expr.dt$top.ind == 0, ]
      top.dt = expr.dt[expr.dt$top.ind == 1, ]
    }
  }
  
  if (model == "Trended") {
    d = DGEList(counts = counts, group = grp.ids)
    d = calcNormFactors(d, method = "RLE")
    e.trd = estimateGLMTrendedDisp(d, design = x)
    trd.fit = glmFit(d, design = x, dispersion = e.trd$trended.dispersion)
    trd.lrt = glmLRT(trd.fit, coef = 2)
    expr.dt$logFC = trd.lrt$table[ ,"logFC"]
    expr.dt$logCPM = trd.lrt$table[ ,"logCPM"]
    expr.dt$FDR = p.adjust(p=trd.lrt$table$PValue, method="BH")
    expr.dt = expr.dt[complete.cases(expr.dt), ]
    # if ask for DE genes by FDR (to be highlighted)
    if (!is.null(FDR)){
      expr.dt$sig.ind = ifelse(test=expr.dt$FDR < FDR, 1, 0)
      colnames(expr.dt) = c("first.grp.mean", "second.grp.mean", "logFC", "logCPM", "FDR", "sig.ind")  
      non.de.dt = expr.dt[expr.dt$sig.ind == 0, ]
      de.dt = expr.dt[expr.dt$sig.ind == 1, ]
    }
    # if ask for top DE genes (to be highlighted)
    if (!is.null(top)){
      expr.dt$top.ind = numeric(dim(expr.dt)[1])
      pv.idx = order(trd.lrt$table$PValue)
      expr.dt = expr.dt[pv.idx, ]
      expr.dt[(1:top),"top.ind"] = 1
      non.top.dt = expr.dt[expr.dt$top.ind == 0, ]
      top.dt = expr.dt[expr.dt$top.ind == 1, ]
    }
  }
  
  ## ------------------------------------------ begin MA plot ------------------------------------------ ##
  
  ## ---------------------------- ##
  ## use FDR cut-off and gradient ##
  ## ---------------------------- ##
  if (gradient & !is.null(FDR)){
    # MA plot with FDR gradients in NBPSeq package
    if ( (model %in% c("NBP","NBQ","NB2")) & !is.null(FDR) ) {
      maplot = 
        ggplot(data=expr.dt, aes(x=logCPM, y=logFC)) + 
        # plot non-significant genes
        geom_point(data=non.de.dt, aes(x=logCPM, y=logFC), colour = "Gray", alpha= 0.9) +
        # plot significant genes (by FDR)
        geom_point(data=de.dt, aes(x=logCPM, y=logFC, colour = FDR), alpha=1) +
        # scaling color gradient for significant genes ONLY
        scale_colour_gradient(name="FDR", limits=c(0,FDR), low = "red", high = "blue") +
        scale_x_log10() +
        #
        geom_hline(yintercept=0, colour="Black", linetype="dashed") + 
        geom_hline(yintercept=-1, colour="Black", linetype="dashed") + 
        geom_hline(yintercept=1, colour="Black", linetype="dashed") + 
        #
        ggtitle(title) + theme_bw() +
        xlab( "Average Expression in RPM (two-group pooled)" ) + 
        ylab(expression(log[2]("Fold Change"))) +
        theme(plot.title = element_text(face = "bold", size = 16), 
              axis.title.x = element_text(face="bold", size=14),
              axis.title.y = element_text(face="bold", size=14, angle=90),
              axis.text.x  = element_text(size=14),
              axis.text.y  = element_text(size=14),
              legend.justification = c(1, 0), 
              legend.position = c(1, 0), 
              legend.title = element_text(face="bold", size=14), 
              legend.key = element_blank(), 
              legend.text = element_text(size = 14))
    }
    
    # MA plot with FDR gradients in edgeR package
    if ( !(model %in% c("NBP","NBQ","NB2")) & !is.null(FDR) ){
      maplot = 
        ggplot(data=expr.dt, aes(x=logCPM, y=logFC)) + 
        # plot non-significant genes
        geom_point(data=non.de.dt, aes(x=logCPM, y=logFC), colour = "Gray", alpha= 0.9) +
        # plot significant genes (by FDR)
        geom_point(data=de.dt, aes(x=logCPM, y=logFC, colour = FDR), alpha=1) +
        # scaling color gradient for significant genes ONLY
        scale_colour_gradient(name="FDR", limits=c(0,FDR), low = "red", high = "blue") +
        #
        geom_hline(yintercept=0, colour="Black", linetype="dashed") + 
        geom_hline(yintercept=-1, colour="Black", linetype="dashed") + 
        geom_hline(yintercept=1, colour="Black", linetype="dashed") + 
        #
        ggtitle(title) + theme_bw() +
        xlab("Average Log Counts Per Million (CPM)") + 
        ylab(expression(log[2]("Fold Change"))) +
        theme(plot.title = element_text(face = "bold", size = 16), 
              axis.title.x = element_text(face="bold", size=14),
              axis.title.y = element_text(face="bold", size=14, angle=90),
              axis.text.x  = element_text(size=14),
              axis.text.y  = element_text(size=14),
              legend.justification = c(1, 0), 
              legend.position = c(1, 0), 
              legend.title = element_text(face="bold", size=14), 
              legend.key = element_blank(), 
              legend.text = element_text(size = 14))
    }
  }
  
  ## ------------------------ ##
  ## use FDR without gradient ##
  ## ------------------------ ##
  if (!gradient & !is.null(FDR)){
    # MA plot with FDR gradients in NBPSeq package
    if ( (model %in% c("NBP","NBQ","NB2")) & !is.null(FDR) ) {
      maplot = 
        ggplot(data=expr.dt, aes(x=logCPM, y=logFC)) + 
        # plot non-significant genes
        geom_point(data=non.de.dt, aes(x=logCPM, y=logFC), colour = "Gray", alpha= 0.9) +
        # plot significant genes (by FDR)
        geom_point(data=de.dt, aes(x=logCPM, y=logFC), colour="Red", alpha=1) +
        scale_x_log10() +
        #
        geom_hline(yintercept=0, colour="Black", linetype="dashed") + 
        geom_hline(yintercept=-1, colour="Black", linetype="dashed") + 
        geom_hline(yintercept=1, colour="Black", linetype="dashed") + 
        #
        ggtitle(title) + theme_bw() +
        xlab( "Average Expression in RPM (two-group pooled)" ) + 
        ylab(expression(log[2]("Fold Change"))) +
        theme(plot.title = element_text(face = "bold", size = 16),
              axis.title.x = element_text(face="bold", size=14),
              axis.title.y = element_text(face="bold", size=14, angle=90),
              axis.text.x  = element_text(size=14),
              axis.text.y  = element_text(size=14),
              legend.position = "none")
    }
    
    # MA plot with FDR gradients in edgeR package
    if ( !(model %in% c("NBP","NBQ","NB2")) & !is.null(FDR) ){
      maplot = 
        ggplot(data=expr.dt, aes(x=logCPM, y=logFC)) + 
        # plot non-significant genes
        geom_point(data=non.de.dt, aes(x=logCPM, y=logFC), colour = "Gray", alpha= 0.9) +
        # plot significant genes (by FDR)
        geom_point(data=de.dt, aes(x=logCPM, y=logFC), colour="Red", alpha=1) +
        #
        geom_hline(yintercept=0, colour="Black", linetype="dashed") + 
        geom_hline(yintercept=-1, colour="Black", linetype="dashed") + 
        geom_hline(yintercept=1, colour="Black", linetype="dashed") + 
        #
        ggtitle(title) + theme_bw() +
        xlab("Average Log Counts Per Million (CPM)") + 
        ylab(expression(log[2]("Fold Change"))) +
        theme(plot.title = element_text(face = "bold", size = 16),
              axis.title.x = element_text(face="bold", size=14),
              axis.title.y = element_text(face="bold", size=14, angle=90),
              axis.text.x  = element_text(size=14),
              axis.text.y  = element_text(size=14),
              legend.position = "none")
    }
  }
  
  # MA plot with TOP genes in NBPSeq package
  if ( (model %in% c("NBP","NBQ","NB2")) & !is.null(top) ) {
    maplot = 
      ggplot(data=expr.dt, aes(x=logCPM, y=logFC)) + 
      # plot non-significant genes
      geom_point(data=non.top.dt, aes(x=logCPM, y=logFC), colour = "Gray", alpha= 0.9) +
      # plot significant genes (by top)
      geom_point(data=top.dt, aes(x=logCPM, y=logFC), colour = "Red", alpha=1) +
      scale_x_log10() +
      #
      geom_hline(yintercept=0, colour="Black", linetype="dashed") + 
      geom_hline(yintercept=-1, colour="Black", linetype="dashed") + 
      geom_hline(yintercept=1, colour="Black", linetype="dashed") + 
      #
      ggtitle(title) + theme_bw() +
      xlab( "Average Expression in RPM (two-group pooled)" ) + 
      ylab(expression(log[2]("Fold Change"))) +
      theme(plot.title = element_text(face = "bold", size = 16), 
            axis.title.x = element_text(face="bold", size=14),
            axis.title.y = element_text(face="bold", size=14, angle=90),
            axis.text.x  = element_text(size=14),
            axis.text.y  = element_text(size=14),
            legend.justification = c(1, 0), 
            legend.position = c(1, 0), 
            legend.title = element_text(face="bold", size=14), 
            legend.key = element_blank(), 
            legend.text = element_text(size = 14))
  }
  
  # MA plot with TOP genes in edgeR package (non-genewise approaches)
  if ( !(model %in% c("NBP","NBQ","NB2","Genewise")) & !is.null(top) ){
    maplot = 
      ggplot(data=expr.dt, aes(x=logCPM, y=logFC)) + 
      # plot non-significant genes
      geom_point(data=non.top.dt, aes(x=logCPM, y=logFC), colour = "Gray", alpha= 0.9) +
      # plot significant genes (by top)
      geom_point(data=top.dt, aes(x=logCPM, y=logFC), colour="Red", alpha=1) +
      #
      geom_hline(yintercept=0, colour="Black", linetype="dashed") + 
      geom_hline(yintercept=-1, colour="Black", linetype="dashed") + 
      geom_hline(yintercept=1, colour="Black", linetype="dashed") + 
      #
      ggtitle(title) + theme_bw() +
      xlab("Average Log Counts Per Million (CPM)") + 
      ylab(expression(log[2]("Fold Change"))) +
      theme(plot.title = element_text(face = "bold", size = 16), 
            axis.title.x = element_text(face="bold", size=14),
            axis.title.y = element_text(face="bold", size=14, angle=90),
            axis.text.x  = element_text(size=14),
            axis.text.y  = element_text(size=14),
            legend.justification = c(1, 0), 
            legend.position = c(1, 0), 
            legend.title = element_text(face="bold", size=14), 
            legend.key = element_blank(), 
            legend.text = element_text(size = 14))
  }
  
  # MA plot with TOP genes in edgeR package (genewise approach only, with dispersion gradient)
  if (gradient & !is.null(top) & model == "Genewise") {
    maplot = 
      ggplot(data=expr.dt, aes(x=logCPM, y=logFC)) + 
      # plot non-significant genes
      geom_point(data=non.top.dt, aes(x=logCPM, y=logFC), colour = "Gray", alpha= 0.9) +
      # plot significant genes (by top)
      geom_point(data=top.dt, aes(x=logCPM, y=logFC, colour = est.disp), alpha=1) +
      # scaling color gradient for significant genes ONLY (varying with estimated dispersions)
      scale_colour_gradient(name="Est.Disp", low = "Blue", high = "Red") +
      #
      geom_hline(yintercept=0, colour="Black", linetype="dashed") + 
      geom_hline(yintercept=-1, colour="Black", linetype="dashed") + 
      geom_hline(yintercept=1, colour="Black", linetype="dashed") + 
      #
      ggtitle(title) + theme_bw() +
      xlab("Average Log Counts Per Million (CPM)") + 
      ylab(expression(log[2]("Fold Change"))) +
      theme(plot.title = element_text(face = "bold", size = 16), 
            axis.title.x = element_text(face="bold", size=14),
            axis.title.y = element_text(face="bold", size=14, angle=90),
            axis.text.x  = element_text(size=14),
            axis.text.y  = element_text(size=14),
            legend.justification = c(1, 0), 
            legend.position = c(1, 0), 
            legend.title = element_text(face="bold", size=14), 
            legend.key = element_blank(), 
            legend.text = element_text(size = 14))
  }
  
  # MA plot with TOP genes in edgeR package (genewise approach only, without dispersion gradient)
  if (!gradient & !is.null(top) & model == "Genewise") {
    maplot = 
      ggplot(data=expr.dt, aes(x=logCPM, y=logFC)) + 
      # plot non-significant genes
      geom_point(data=non.top.dt, aes(x=logCPM, y=logFC), colour = "Gray", alpha= 0.9) +
      # plot significant genes (by top)
      geom_point(data=top.dt, aes(x=logCPM, y=logFC), colour="Red", alpha=1) +
      #
      geom_hline(yintercept=0, colour="Black", linetype="dashed") + 
      geom_hline(yintercept=-1, colour="Black", linetype="dashed") + 
      geom_hline(yintercept=1, colour="Black", linetype="dashed") + 
      #
      ggtitle(title) + theme_bw() +
      xlab("Average Log Counts Per Million (CPM)") + 
      ylab(expression(log[2]("Fold Change"))) +
      theme(plot.title = element_text(face = "bold", size = 16),
            axis.title.x = element_text(face="bold", size=14),
            axis.title.y = element_text(face="bold", size=14, angle=90),
            axis.text.x  = element_text(size=14),
            axis.text.y  = element_text(size=14),
            legend.position = "none")
  }
  
  return(maplot)
  
}
