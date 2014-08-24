

#' @title Fitting dispersion models and getting p-values for the power-robustness simulations
#' 
#' @details We use the normalization method discussed in Anders and Huber (2010) for all dispersion models.
#' For QL and its variants, we use the upper 75% quantile approach as suggested in the package manual.
#' 
#' @param sim.data.obj  output from sim.data.prepare()
#' 
#' @param model  NB dispersion model to be fitted, including: "Common", "NBP", "NBQ", "NBS", "Trended", "Genewise", "Tagwise-Common", 
#' "Tagwise-Trend", "QL", "QLShrink", "QLSpline" and "Genewise-HOA"
#' 
#' @return a list of two elements: 
#' (1) dispersion model name ("model"); 
#' (2) DE test p-value, used for calculating TPR-FDR curve  ("pval.y")
#' 
#' @export
#' 
#' @author Gu Mi
#' 
sim.model.check = function(sim.data.obj, model=NULL){
  
  y = sim.data.obj$counts
  x = sim.data.obj$design
  grp.ids = factor(apply(x, 1, function(x) { paste(rev(x), collapse = ".") }), labels = seq(ncol(x)))
  
  if (is.null(model)){
    stop("You must specify one of the following dispersion models: 'Common', 'NBP', 'NBQ', 'NBS', 'Trended', 'Genewise', 'Tagwise-Common', 'Tagwise-Trend', 'QL', 'QLShrink', 'QLSpline', or 'Genewise-HOA'!")
  }
  stopifnot(model %in% c("Common","NBP", "NBQ", "NBS", "Trended", "Genewise","Tagwise-Common", "Tagwise-Trend", 'QL', 'QLShrink', 'QLSpline', 'Genewise-HOA'))
  
  if (model == "NBP"){
    nbp.fit = nb.glm.test(counts=y, x=x, beta0=c(NA,0), normalization.method="AH2010", dispersion.model = "NBP", tests="HOA")    
    pval.y = nbp.fit$test.results$HOA$p.values
  }
  
  if (model == "NBQ"){
    nbq.fit = nb.glm.test(counts=y, x=x, beta0=c(NA,0), normalization.method="AH2010", dispersion.model = "NBQ", tests="HOA")
    pval.y = nbq.fit$test.results$HOA$p.values
  }
  
  if (model == "NBS"){
    nbs.fit = nb.glm.test(counts=y, x=x, beta0=c(NA,0), normalization.method="AH2010", dispersion.model = "NBS", tests="HOA")
    pval.y = nbs.fit$test.results$HOA$p.values
  }
  
  if (model == "Common"){
    d = DGEList(counts = y, group = grp.ids)
    # Normalize: method="RLE" is the scaling factor method proposed by Anders and Huber (2010).
    d = calcNormFactors(d, method="RLE")
    e.com = estimateGLMCommonDisp(d, x)   
    com.fit = glmFit(d, x, dispersion=e.com$common.dispersion)
    com.lrt = glmLRT(com.fit,coef=2)
    pval.y = com.lrt$table$PValue
  }
  
  # ------------------------------------------------------------------------- #
  # Genewise-HOA test 
  # ------------------------------------------------------------------------- #
  
  if (model == "Genewise-HOA") {
    nf = estimate.norm.factors(counts=y, method="AH2010")
    nb.data = prepare.nb.data(counts=y, norm.factors=nf)
    grp1 = as.character(unique(grp.ids)[1])
    grp2 = as.character(unique(grp.ids)[2])
    res = genewise.hoa(nb.data=nb.data, grp.ids=grp.ids, grp1=grp1, grp2=grp2, R = 100, print.level=1)
    
    #names(res)
    
#     DEBUG = FALSE
#     if (DEBUG){
#       library(SeqDisp)
#       data(mouse)
#       dt = mouse
#       cpms = cpm(dt)
#       keep = rowSums(cpms > 1) >= 1
#       dt2 = dt[keep, ]
#       nf = estimate.norm.factors(dt2, method = "AH2010")
#       set.seed(539)
#       idx = sample(1:dim(dt2)[1], size = 5000, replace = FALSE)
#       y = dt2[idx,1:6]
#       n=6
#       grp.ids = as.factor(c(rep(1, n/2), rep(2, n/2)))      
#     }
    
    # Pay attention to the NA p-values for HOA tests: now assign 1 for those NA p-values
    pval.y = res$p.hoa
    pval.y[is.na(pval.y)] = 1
  }  
  
#   if (model == "Genewise"){
#     d = DGEList(counts = y, group = grp.ids)
#     d = calcNormFactors(d, method="RLE")
#     e.com = estimateGLMCommonDisp(d, design=x, verbose=FALSE) 
#     e.gen = estimateGLMTagwiseDisp(d, design=x, dispersion=e.com$common.dispersion, prior.df=0, trend=FALSE)  # prior.df = 0
#     # e.gen: trend is FALSE since we didn't use estimateGLMTrendedDisp() beforehand
#     gen.fit = glmFit(d, design=x, dispersion=e.gen$tagwise.dispersion)
#     gen.lrt = glmLRT(gen.fit,coef=2)
#     pval.y = gen.lrt$table$PValue
#   }
  
  if (model == "Tagwise-Common"){
    d = DGEList(counts = y, group = grp.ids)
    d = calcNormFactors(d, method="RLE")
    d = estimateGLMCommonDisp(d, x)   
    e.tgc = estimateGLMTagwiseDisp(d, design=x, trend=FALSE)
    tgc.fit = glmFit(d, x, dispersion=e.tgc$tagwise.dispersion)
    tgc.lrt = glmLRT(tgc.fit,coef=2)
    pval.y = tgc.lrt$table$PValue
  }
  
  if (model == "Tagwise-Trend"){
    d = DGEList(counts = y, group = grp.ids)
    d = calcNormFactors(d, method="RLE")
    d = estimateGLMTrendedDisp(d, design=x)
    e.tgt = estimateGLMTagwiseDisp(d, design=x, trend=TRUE)
    tgt.fit = glmFit(d, design=x, dispersion=e.tgt$tagwise.dispersion)
    tgt.lrt = glmLRT(tgt.fit,coef=2)
    pval.y = tgt.lrt$table$PValue
  }
  
  if (model == "Trended"){
    d = DGEList(counts = y, group = grp.ids)
    d = calcNormFactors(d, method="RLE")
    e.trd = estimateGLMTrendedDisp(d, design=x)
    trd.fit = glmFit(d, design=x, dispersion=e.trd$trended.dispersion)
    trd.lrt = glmLRT(trd.fit,coef=2)
    pval.y = trd.lrt$table$PValue
  }
  
  if (model == "QL"){
    design.list = vector("list",2)
    design.list[[1]] = x 
    design.list[[2]] = rep(1,length(grp.ids))
    log.offset = log(apply(y,2,quantile,0.75))
    fit = QL.fit(counts=y, design.list=design.list, log.offset=log.offset, Model="NegBin",method="glm", NBdisp="trend", print.progress=FALSE)
    res = QL.results(fit, Plot=FALSE)
    pval.y = as.numeric(res$P.values$QL)
  }
  
  if (model == "QLShrink"){
    design.list = vector("list",2)
    design.list[[1]] = x  
    design.list[[2]] = rep(1,length(grp.ids))
    log.offset = log(apply(y,2,quantile,0.75))
    fit = QL.fit(counts=y, design.list=design.list, log.offset=log.offset, Model="NegBin",method="glm", NBdisp="trend", print.progress=FALSE)
    res = QL.results(fit, Plot=FALSE)
    pval.y = as.numeric(res$P.values$QLShrink)
  }
  
  if (model == "QLSpline"){
    design.list = vector("list",2)
    design.list[[1]] = x  
    design.list[[2]] = rep(1,length(grp.ids))
    log.offset = log(apply(y,2,quantile,0.75))
    fit = QL.fit(counts=y, design.list=design.list, log.offset=log.offset, Model="NegBin",method="glm", NBdisp="trend", print.progress=FALSE)
    res = QL.results(fit, Plot=FALSE)
    pval.y = as.numeric(res$P.values$QLSpline)
  }
  
  return(list(model=model, pval.y=pval.y))
}

