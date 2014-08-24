
#' @title Preparing data frames for producing mean-dispersion plot using ggplot2 for several NB dispersion models
#' 
#' @details We temporarily use method-of-moments (MOM) estimators for each gene in the scatter plot.
#' 
#' @param counts  RNA-Seq raw read count matrix
#' 
#' @param x  design matrix
#' 
#' @param model  name of the NB dispersion model, including "NBP", "NBQ", "NBS", "STEP", "Common", 
#' "Tagwise-Common", "Tagwise-Trend" and "Trended"
#' 
#' @param sc  a scaling factor to MOM estimator phi, for an overall up-/down-shift of the dispersions (default = 1)
#' 
#' @return depending on the model, this function returns coordinates of the scatter plot and names of the dispersion models
#' 
#' @author Gu Mi
#' 
#' @keywords internal
#' 
mddata = function(counts, x, model = NULL, sc = 1){
  
  m = dim(counts)[1]
  n = dim(counts)[2]
  
  # naive estimations of relative frequency and phi for the basic scatter plot
  # it makes sense to have only one group (design matrix x is a column: intercept-only model)
  
  # method of moments (MOM) estimates for mu
  mu.mom = edgeR:::expandAsMatrix(rowMeans(counts), dim=c(m,n))
  
  # relative frequencies pi = mu / N
  re.freq = (mu.mom / (matrix(1, m, 1) %*% matrix(colSums(counts), 1, n)))[ ,1] 
  #qt.re.freq = quantile(re.freq, c(0.001, 0.999))
  phi.mom = sc * ( rowSums((counts - mu.mom)^2) - rowSums(mu.mom) ) / rowSums(mu.mom^2)  # may have NaN
  #id = (phi.mom > 0 & !is.nan(phi.mom) & re.freq > qt.re.freq[1] & re.freq < qt.re.freq[2])  
  id = (phi.mom > 0 & !is.nan(phi.mom))
  # may discard some phi.hat here not in the plotting: this "id" is used "globally" to subset
  n.pt = sum(id)
  
  if (is.null(model)){
    coords.scatter = data.frame(x0=re.freq[id], y0=phi.mom[id])   # for scatter plot
    return(coords.scatter)
  }
  
  if (model == "NBP"){
    nf = estimate.norm.factors(counts=counts, lib.sizes=colSums(counts), method="AH2010")
    nb.data = prepare.nb.data(counts, lib.sizes=colSums(counts), norm.factors=nf)
    nbp.disp.z = NBPSeq:::disp.nbp(counts=counts, eff.lib.sizes=nb.data$eff.lib.sizes, x=x)   # to get "z"
    nbp.disp = estimate.dispersion(nb.data = nb.data, x = x, model = "NBP", method = "MAPL")    # to get disp. estimates
    phi.nbp = nbp.disp$estimates
    pi.hat = nbp.disp.z$pi.pre[,1]
    coords = data.frame(x=pi.hat, y=phi.nbp[ ,1])
    return(as.data.frame(cbind(Dispersion.Model = rep("NBP", m), coords)))
  }
  
  if (model == "NBQ"){
    nf = estimate.norm.factors(counts=counts, lib.sizes=colSums(counts), method="AH2010")
    nb.data = prepare.nb.data(counts, lib.sizes=colSums(counts), norm.factors=nf)
    nbq.disp.z = NBPSeq:::disp.nbq(counts=counts, eff.lib.sizes=nb.data$eff.lib.sizes, x=x)   # to get "z"
    nbq.disp = estimate.dispersion(nb.data = nb.data, x = x, model = "NBQ", method = "MAPL")  # to get disp. estimates
    phi.nbq = nbq.disp$estimates
    pi.hat = nbq.disp.z$pi.pre[,1]
    coords = data.frame(x=pi.hat, y=phi.nbq[ ,1])
    return(as.data.frame(cbind(Dispersion.Model = rep("NBQ", m), coords)))  
  }
  
  if (model == "NBS"){
    nf = estimate.norm.factors(counts=counts, lib.sizes=colSums(counts), method="AH2010")
    nb.data = prepare.nb.data(counts, lib.sizes=colSums(counts), norm.factors=nf)
    nbs.disp.z = NBPSeq:::disp.nbs(counts=counts, eff.lib.sizes=nb.data$eff.lib.sizes, x=x)   # to get "z"
    nbs.disp = estimate.dispersion(nb.data = nb.data, x = x, model = "NBS", method = "MAPL")
    phi.nbs = nbs.disp$estimates
    pi.hat = nbs.disp.z$pi.pre[,1]
    coords = data.frame(x=pi.hat, y=phi.nbs[ ,1])
    return(as.data.frame(cbind(Dispersion.Model = rep("NBS", m), coords)))  
  }
  
  if (model == "STEP"){
    nf = estimate.norm.factors(counts=counts, lib.sizes=colSums(counts), method="AH2010")
    nb.data = prepare.nb.data(counts, lib.sizes=colSums(counts), norm.factors=nf)
    nbstep.disp.z = NBPSeq:::disp.step(counts=counts, eff.lib.sizes=nb.data$eff.lib.sizes, x=x, df=10)   # to get "z"
    nbstep.disp = estimate.dispersion(nb.data = nb.data, x = x, model = "NB2", method = "MAPL")
    phi.nbstep = nbstep.disp$estimates
    pi.hat = nbstep.disp.z$pi.pre[,1]
    coords = data.frame(x=pi.hat, y=phi.nbstep[ ,1])
    return(as.data.frame(cbind(Dispersion.Model = rep("STEP", m), coords)))  
  }
  
  if (model == "Common"){
    y.dge = DGEList(counts)
    y.dge = calcNormFactors(y.dge, method = "RLE")
    e.com = estimateGLMCommonDisp(y.dge, x)
    phi = rep(e.com$common.dispersion, m)
    rel.freq = exp(e.com$AveLogCPM*log(2) - log(1e+06))
    # for plotting purpose:
    id2 = order(rel.freq)
    rel.freq.ord = rel.freq[id2]
    phi.ord = phi[id2]
    coords = data.frame(x=rel.freq.ord, y=phi.ord)   
    return(as.data.frame(cbind(Dispersion.Model = rep("Common", m), coords)))
  }
  
  if (model == "Tagwise-Common"){
    y.dge = DGEList(counts)
    y.dge = calcNormFactors(y.dge, method = "RLE")
    y.dge = estimateGLMCommonDisp(y.dge, x)
    e.tgc = estimateGLMTagwiseDisp(y.dge, design=x, trend=FALSE)
    phi = e.tgc$tagwise.dispersion
    rel.freq = exp(e.tgc$AveLogCPM*log(2) - log(1e+06))
    # for plotting purpose:
    id2 = order(rel.freq)
    rel.freq.ord = rel.freq[id2]
    phi.ord = phi[id2]
    coords = data.frame(x=rel.freq.ord, y=phi.ord)   
    return(as.data.frame(cbind(Dispersion.Model = rep("Tagwise-Common", m), coords)))
  }
  
  if (model == "Tagwise-Trend"){
    y.dge = DGEList(counts)
    y.dge = calcNormFactors(y.dge, method = "RLE")
    y.dge = estimateGLMTrendedDisp(y.dge, x)
    e.tgt = estimateGLMTagwiseDisp(y.dge, design=x, trend=TRUE)
    phi = e.tgt$tagwise.dispersion
    rel.freq = exp(e.tgt$AveLogCPM*log(2) - log(1e+06))
    # for plotting purpose:
    id2 = order(rel.freq)
    rel.freq.ord = rel.freq[id2]
    phi.ord = phi[id2]
    coords = data.frame(x=rel.freq.ord, y=phi.ord)    
    return(as.data.frame(cbind(Dispersion.Model = rep("Tagwise-Trend", m), coords))) 
  }
  
  if (model == "Trended"){ 
    trd = dispBinTrend(counts)
    # names(trd)  # "AveLogCPM"      "dispersion"     "bin.AveLogCPM"  "bin.dispersion"
    phi = trd$dispersion
    rel.freq = exp(trd$AveLogCPM*log(2) - log(1e+06))
    # for plotting purpose:
    id2 = order(rel.freq)
    rel.freq.ord = rel.freq[id2]
    phi.ord = phi[id2]
    coords = data.frame(x=rel.freq.ord, y=phi.ord)
    return(as.data.frame(cbind(Dispersion.Model = rep("Trended", m), coords)))
  }
}




#' @title Plot mean-dispersion scatter plot (log-log scale) for several NB dispersion models using ggplot2
#' 
#' @details We temporarily use method-of-moments (MOM) estimators for each gene in the scatter plot, and superimpose
#' fitted curves for the dispersion models specified in \code{model.vec}.
#' 
#' @param model.vec  a vector of NB dispersion models where the fitted curves will be superimposed on the scatter plot. 
#' Specify: "NBP", "NBQ", "NBS", "STEP", "Common", "Tagwise-Common", "Tagwise-Trend" or "Trended"
#' 
#' @param counts  RNA-Seq raw read count matrix
#' 
#' @param x  design matrix
#' 
#' @param sc  a scaling factor to MOM estimator phi, for an overall up-/down-shift of the dispersions (default = 1)
#' 
#' @param title  main title of the plot
#' 
#' @param data.note  note of the dataset (put in the lower-left corner of the plot)
#' 
#' @return depending on the model, this function returns coordinates of the scatter plot and names of the dispersion models
#' 
#' @author Gu Mi
#' 
#' @export
#' 
pl.md = function(model.vec, counts, x, sc = 1, title=NULL, data.note=NULL){
  
  k = length(model.vec)
  result.lst = rep( list(NA), k) 
  for (i in seq_len(k)){
    result.lst[[i]] = mddata(counts, x, model = model.vec[i])
  }
  df.long = do.call(rbind.data.frame, result.lst)
  n.models = length(model.vec)
  
  dta = mddata(counts, x, model=NULL, sc = sc)
  
  md = 
    ggplot(data = dta, aes(x = x0, y = y0)) + 
    # scatter plot of MOM estimates
    geom_point(data = dta, aes(x = x0, y = y0), alpha=I(0.25), size=I(1)) +    
    # x axis
    scale_x_log10("Estimated Relative Frequency", breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)), limits=c(min(dta$x0), max(dta$x0))) + 
    # y axis
    scale_y_log10("Estimated NB Dispersion Parameter", breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)), limits=c(min(dta$y0), max(dta$y0))) +
    # fitted curves
    geom_line(data = df.long, aes(x=x, y=y, colour=factor(Dispersion.Model), size=factor(Dispersion.Model),
                                  linetype=factor(Dispersion.Model), alpha=factor(Dispersion.Model))) +
    scale_size_manual(values = rep(1.2, n.models)) + 
    scale_linetype_manual(values = seq(1, length(unique(df.long$Dispersion.Model)))) + 
    scale_alpha_manual(values = rep(0.8, n.models)) + 
    #
    theme_bw() +
    ggtitle(title) + 
    # put data.note at lower-left corner
    annotate("text", x=min(dta$x0), y=min(dta$y0), hjust=0, vjust=0, label=data.note, size = 5, fontface = 3) +
    theme(plot.title = element_text(face="bold", size=16),
          axis.title.x = element_text(face="bold", size=16),
          axis.title.y = element_text(face="bold", size=16, angle=90),
          axis.text.x  = element_text(size=12),
          axis.text.y  = element_text(size=12),
          legend.key.width = unit(3, "line"),
          legend.justification=c(1,0), 
          legend.position=c(1,0),
          legend.title = element_blank(),
          legend.key = element_blank(),
          legend.text = element_text(size=14)
          )
}