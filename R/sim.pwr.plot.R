
#' @title Preparing and plotting all dispersion models using ggplot2
#' 
#' @details This function prepares meta-data and calls the following functions to obtain ggplot2 results:
#'    
#'    tpr.fdr.plot(): for TPR-FDR plot
#'    n.tp.plot(): for true positive plot (among top n genes)
#'    prec.plot(): for Precision-Recall plot 
#' 
#' @param sim.data.obj  output from sim.data.prepare().
#' 
#' @param type  which type of plot to draw? Either "TPR-FDR", "N-TP" or "PREC"
#' 
#' @param beta  a non-negative real value for calculating the general F_beta measure (default = 1: recall and precision are evenly weighted)
#' 
#' @param model.vec  a character vector specifying the names of the NB dispersion models. Currently supported include 
#' "NBP", "NBQ", "NBS", "Common", "Genewise", "Tagwise-Common", "Tagwise-Trend", "Trended", "QL", "QLShrink", "QLSpline" and "Genewise-HOA".
#' 
#' @param title.note  extra data annotation in parentheses.
#' 
#' @param legend  whether to put a legend on a single plot.
#' 
#' @param legend.just  the legend.justification argument in ggplot2 controling the position of the legend.
#'
#' @param N  the number of top genes selected
#'
#' @param x.limit  specify the plotting range for the x-axis (TPR), for better visualization
#' 
#' @param panel.label  the label on the lower-right corner of each plot (for journal figure legends).
#' 
#' @return a ggplot object of the specified curves. Use print() to draw the specified curves.
#' 
#' @export
#' 
#' @author Gu Mi
#' 
sim.pwr.plot = function(sim.data.obj, type="TPR-FDR", beta=1, model.vec=NULL, title.note="", legend=TRUE, 
                        legend.just=c(0,1), N=500, x.limit=c(0,1), panel.label="A"){
  
  m = sim.data.obj$m
  idx.DE = sim.data.obj$idx.DE
  true.de.vec = logical(length=m)  # all FALSE
  true.de.vec[idx.DE] = TRUE       # the "grp" argument in pwrdata(): this is the true DE label vector
  
  model.lst = rep( list(NA), length(model.vec)) 
  for (i in seq_len(length(model.vec))){
    model.lst[[i]] = sim.model.check(sim.data.obj=sim.data.obj, model=model.vec[i])
  }
  
  n.models = length(model.lst)  # number of models compared
  all_methods = rep( list(NA), n.models) 
  model.names = rep( "", n.models)
  for (i in 1:n.models){
    all_methods[[i]] = data.frame(grp = true.de.vec, pvals = model.lst[[i]]$pval.y)   # grp (true DE labels); pvals (DE p-values)
    model.names[[i]] = model.lst[[i]]$model
  }
  names(all_methods) = model.names
  
  # Each data frame containing our data is then put together into a list object which we pass to the tpr.fdr.plot() or prec.plot()
  #
  if (type == "TPR-FDR") {
    p = tpr.fdr.plot(all_methods, title = title.note, legend=legend, legend.just=legend.just, x.limit=x.limit, panel.label=panel.label)
    return(p)  # a list of two components, pwr and info
  }
  
  if (type == "N-TP"){
    p = n.tp.plot(all_methods, title = title.note, legend=legend, legend.just=legend.just, N=N, panel.label=panel.label)
    return(p)  # a list of two components, pwr and info
  }
  
  if (type == "PREC") {
    p = prec.plot(all_methods, beta=beta, title = title.note, legend=legend, legend.just=legend.just, x.limit=x.limit, panel.label=panel.label)
    return(p)  # a list of three components, pwr, info and Fmeasure
  }
  
}


##' @title Preparing quantities for TPR vs. FDR plot, or Precision-Recall plot
##' 
##' @details This function produces x and y co-ordinates for either type of plots.
##' 
##' @param grp  labels classifying subject status (true DE labels)
##' 
##' @param pred  values of each observation (DE test p-values)
##' 
##' @param type  which type of plot to draw? Either "TPR-FDR", "N-TP" or "PREC"
##' 
##' @param N  the number of top genes selected
##' 
##' @param beta  a non-negative real value for calculating the general F_beta measure (default = 1: recall and precision are evenly weighted)
##' 
##' @return List with the common component: 
##' 
##' pwr: data.frame with x and y co-ordinates of plot (TPR and FDR, or recall and precision), 
##'      and sorted p-values, q-values by B-H adjustment and F-measures (if PR-curve is desired).
##' 
##' @keywords internal
##' 
##' @author Gu Mi
##' 
pwrdata = function(grp, pred, type="TPR-FDR", N = 500, beta=beta){ 
  
  grp = as.factor(grp)  # a factor w/ levels "FALSE" "TRUE": 1 1 1 1 1 1 1 1 1 2 ...
  
  if (length(pred) != length(grp)) {
    stop("The number of classifiers must match the number of data points")
  } 
  if (length(levels(grp)) != 2) {
    stop("There must only be 2 values for the classifier")
  }
  
  #cut = unique(pred)  # a set of thresholds
  tmp1 = cbind(grp, pred)          # for grp: 1 = FALSE; 2 = TRUE
  tmp2 = tmp1[order(tmp1[ ,2]), ]  # rank the p-values (pred)
  grp = as.factor(tmp2[ ,1])       # the corresponding grp after ranking the p-values (pred)
  cut = pred = tmp2[ ,2]   # a set of thresholds (p-values) from smallest (~0) to largest (~1)
  
  # true positives: for every cutoff value we count # of entries where the predictor vector (DE p-values) is less than the cutoff AND 
  # the group vector is the higher factor level ("TRUE", that's why we ensure it is an ordered factor)
  # levels(grp)[1] --> "1" (FALSE)
  # levels(grp)[2] --> "2" (TRUE)
  tp = sapply(cut, function(x) length(which(pred < x & grp == levels(grp)[2])))
  fn = sapply(cut, function(x) length(which(pred > x & grp == levels(grp)[2])))
  fp = sapply(cut, function(x) length(which(pred < x & grp == levels(grp)[1])))
  tn = sapply(cut, function(x) length(which(pred > x & grp == levels(grp)[1])))
  q.value = p.adjust(pred, method="BH")
  
  if (type == "TPR-FDR") {
    tpr = tp / (tp + fn)        # True Positive Rate (x-axis) --> power
    fdr = fp / (tp + fp)        # False Discovery Rate (y-axis)
    pwr = data.frame(x = tpr, y = fdr, p = pred, q = q.value)  
    # does NOT matter if sort tpr or not, as power is definitely an increasing function
    # pwr = pwr[order(pwr$x, pwr$y), ]
    pwr = pwr[complete.cases(pwr), ]
    return(list(pwr=pwr))
  }
  
  if (type == "N-TP") {
    pwr = data.frame(x = (1:length(pred))[1:N], y = tp[1:N], p = pred, q = q.value)
    pwr = pwr[complete.cases(pwr), ]
    return(list(pwr=pwr))
  }
  
  if (type == "PREC") {
    recall = tp / (tp + fn)     # True Positive Rate (sensitivity), or Recall (x-axis)
    precision = tp / (tp + fp)  # Precision (y-axis)
    Fmeasure = (1 + beta^2) * precision * recall / ( (beta^2) * precision + recall )
    pwr = data.frame(x = recall, y = precision, p = pred, q = q.value, Fmeas = Fmeasure)
    pwr = pwr[order(pwr$x, pwr$y), ]
    pwr = pwr[complete.cases(pwr), ]
    return(list(pwr=pwr))
  }
  
}




##' @title Plot the True Positve Rate (TPR) vs. False Discovery Rate (FDR) plot using ggplot2 (with other information reported)
##' 
##' @details This function produces a single ggplot2 object, and is made invisible to end-users. For multiple plots, see the
##' \code{pl.tpr.fdr} function. 
##' 
##' @param test.data.list  meta-data of the quantities for each dispersion model, as prepared in sim.pwr.plot()
##' 
##' @param groupName  leave as default "grp"; c.f. sim.pwr.plot()
##' 
##' @param predName  leave as default "pvals"; c.f. sim.pwr.plot()
##' 
#' @param legend  whether to put a legend on a single plot
#' 
#' @param legend.just  the legend.justification argument in ggplot2 controling the position of the legend
#' 
#' @param x.limit  specify the plotting range for the x-axis (TPR), for better visualization
#' 
#' @param panel.label  the label on the lower-right corner of each plot (for journal figure legends)
#' 
#' @param title  the plot title for a single panel
#' 
#' @return a list of a ggplot2 object and relevant information. 
#' 
#' @keywords internal
#' 
#' @importFrom plyr llply
#' @importFrom plyr ldply
#' 
#' @author Gu Mi
#' 
tpr.fdr.plot = function(test.data.list, groupName = "grp", predName = "pvals", legend=TRUE, 
                        legend.just=c(0,1), x.limit=c(0,1), panel.label="A", title = "Plot Title") {
  
  plotdata = llply(test.data.list, function(x) 
    with(x, pwrdata(grp = eval(parse(text = groupName)), pred = eval(parse(text = predName)), type="TPR-FDR") ))
  plotdata2 = list(pwr = ldply(plotdata, function(x) x$pwr))
  names(plotdata2$pwr) = c("model.names", "tpr", "fdr", "pvalue", "qvalue")
  # by default, ".id" is used for the column of model names
  
  if (legend) {
    p = 
      ggplot(plotdata2$pwr, aes(x = tpr, y = fdr)) +
      geom_line(size=0.5, aes(linetype=factor(model.names), colour=factor(model.names))) +
      theme_bw() +
      scale_x_continuous("", limits=x.limit) +
      scale_y_continuous("", limits=c(0,1)) +
      # in order to get legend of lines having both colour and linetype, they must have the same title and labels, etc.
      # o.w. there will be duplicate legends produced
      scale_linetype_discrete(names(test.data.list),
                              breaks = names(test.data.list),
                              labels = names(test.data.list)) +
      scale_colour_brewer(palette="Dark2", 
                          names(test.data.list),
                          breaks = names(test.data.list),
                          labels = names(test.data.list)) +
      ggtitle(title) + 
      annotate("text", x = x.limit[2], y = 0, label=panel.label, size = 10, face="bold") +
      theme(plot.title = element_text(face="bold", size=18), 
            axis.title.x = element_text(face="bold", size=16),
            axis.title.y = element_text(face="bold", size=16, angle=90),
            axis.text.x  = element_text(size=16),
            axis.text.y  = element_text(size=16),
            plot.margin=unit(x=c(0.8, 0.8, 1, 1), units="cm"),  # top, right, bottom, and left
            legend.key.width = unit(3, "line"),
            legend.justification = legend.just, 
            legend.position = legend.just,
            legend.title = element_blank(),
            legend.key = element_blank(),
            legend.text = element_text(size=16)
      )
  }
  #
  if (!legend) {
    p = ggplot(plotdata2$pwr, aes(x = tpr, y = fdr)) +
      geom_line(size=0.5, aes(linetype=factor(model.names), colour=factor(model.names))) +
      theme_bw() +
      scale_x_continuous("", limits=x.limit) +
      scale_y_continuous("", limits=c(0,1)) +
      # in order to get legend of lines having both colour and linetype, they must have the same title and labels, etc.
      # o.w. there will be duplicate legends produced
      scale_linetype_discrete(names(test.data.list),
                              breaks = names(test.data.list),
                              labels = names(test.data.list)) +
      scale_colour_brewer(palette="Dark2",
                          names(test.data.list),
                          breaks = names(test.data.list),
                          labels = names(test.data.list)) +
      ggtitle(title) + 
      annotate("text", x = x.limit[2], y = 0, label=panel.label, size = 10, face="bold") +
      theme(plot.title = element_text(face="bold", size=18), 
            axis.title.x = element_text(face="bold", size=16),
            axis.title.y = element_text(face="bold", size=16, angle=90),
            axis.text.x  = element_text(size=16),
            axis.text.y  = element_text(size=16),
            plot.margin=unit(x=c(0.8, 0.8, 1, 1), units="cm"),  # top, right, bottom, and left
            legend.position = "none"
      )
  }
  return(list(pwr = p, info = plotdata2$pwr))
}


##' @title Plot the top N genes vs. True Positve Rate (TPR) plot using ggplot2 (with other information reported)
##' 
##' @details This function produces a single ggplot2 object, and is made invisible to end-users. For multiple plots, see the
##' \code{pl.n.tp} function. 
##' 
##' @param test.data.list  meta-data of the quantities for each dispersion model, as prepared in sim.pwr.plot()
##' 
##' @param groupName  leave as default "grp"; c.f. sim.pwr.plot()
##' 
##' @param predName  leave as default "pvals"; c.f. sim.pwr.plot()
##' 
#' @param legend  whether to put a legend on a single plot
#' 
#' @param legend.just  the legend.justification argument in ggplot2 controling the position of the legend
#' 
#' @param N  the number of top genes selected
#'  
#' @param panel.label  the label on the lower-right corner of each plot (for journal figure legends)
#' 
#' @param title  the plot title for a single panel
#' 
#' @return a list of a ggplot2 object and relevant information. 
#' 
#' @keywords internal
#' 
#' @importFrom plyr llply
#' @importFrom plyr ldply
#' 
#' @author Gu Mi
#' 
n.tp.plot = function(test.data.list, groupName = "grp", predName = "pvals", legend=TRUE,
                     legend.just=c(0,1), N=500, panel.label="A", title = "Plot Title") {
  
  plotdata = llply(test.data.list, function(x) 
    with(x, pwrdata(grp = eval(parse(text = groupName)), pred = eval(parse(text = predName)), type="N-TP", N=N) ))
  plotdata2 = list(pwr = ldply(plotdata, function(x) x$pwr))
  names(plotdata2$pwr) = c("model.names", "N", "tp", "pvalue", "qvalue")
  # by default, ".id" is used for the column of model names
  
  if (legend) {
    p = 
      ggplot(plotdata2$pwr, aes(x = N, y = tp)) +
      geom_line(size=0.5, aes(linetype=factor(model.names), colour=factor(model.names))) +
      theme_bw() +
      scale_x_continuous("", limits=c(1,N)) +
      scale_y_continuous("", limits=c(0, max(N/2, max(plotdata2$pwr[1:N,"tp"])))) +
      #scale_x_continuous("") +
      #scale_y_continuous("") +
      # in order to get legend of lines having both colour and linetype, they must have the same title and labels, etc.
      # o.w. there will be duplicate legends produced
      scale_linetype_discrete(names(test.data.list),
                              breaks = names(test.data.list),
                              labels = names(test.data.list)) +
      scale_colour_brewer(palette="Dark2", 
                          names(test.data.list),
                          breaks = names(test.data.list),
                          labels = names(test.data.list)) +
      ggtitle(title) + 
      annotate("text", x = N, y = 1, label=panel.label, size = 9, face="bold") +
      theme(plot.title = element_text(face="bold", size=16), 
            axis.title.x = element_text(face="bold", size=14),
            axis.title.y = element_text(face="bold", size=14, angle=90),
            axis.text.x  = element_text(size=14),
            axis.text.y  = element_text(size=14),
            plot.margin=unit(x=c(0.8, 0.8, 1, 1), units="cm"),  # top, right, bottom, and left
            legend.key.width = unit(3, "line"),
            legend.justification = legend.just, 
            legend.position = legend.just,
            legend.title = element_blank(),
            legend.key = element_blank(),
            legend.text = element_text(size=14)
      )
  }
  #
  if (!legend) {
    p = 
      ggplot(plotdata2$pwr, aes(x = N, y = tp)) +
      geom_line(size=0.5, aes(linetype=factor(model.names), colour=factor(model.names))) +
      theme_bw() +
      scale_x_continuous("", limits=c(1,N)) +
      scale_y_continuous("", limits=c(0, max(N/2, max(plotdata2$pwr[1:N,"tp"])))) +
      #scale_x_continuous("") +
      #scale_y_continuous("") +
      # in order to get legend of lines having both colour and linetype, they must have the same title and labels, etc.
      # o.w. there will be duplicate legends produced
      scale_linetype_discrete(names(test.data.list),
                              breaks = names(test.data.list),
                              labels = names(test.data.list)) +
      scale_colour_brewer(palette="Dark2",
                          names(test.data.list),
                          breaks = names(test.data.list),
                          labels = names(test.data.list)) +
      ggtitle(title) + 
      annotate("text", x = N, y = 1, label=panel.label, size = 9, face="bold") +
      theme(plot.title = element_text(face="bold", size=16), 
            axis.title.x = element_text(face="bold", size=14),
            axis.title.y = element_text(face="bold", size=14, angle=90),
            axis.text.x  = element_text(size=14),
            axis.text.y  = element_text(size=14),
            plot.margin=unit(x=c(0.8, 0.8, 1, 1), units="cm"),  # top, right, bottom, and left
            legend.position = "none"
      )
  }
  return(list(pwr = p, info = plotdata2$pwr))
}



##' @title Plot the Precision-Recall plot using ggplot2 (with other information and F-measures reported)
##' 
##' @details This function produces a single ggplot2 object, and is made invisible to end-users. For multiple plots, see the
##' \code{pl.prec} function. 
##' 
##' @param test.data.list  meta-data of the quantities for each dispersion model, as prepared in sim.pwr.plot()
##' 
##' @param beta  a non-negative real value for calculating the general F_beta measure (default = 1: recall and precision are evenly weighted)
##' 
##' @param groupName  leave as default "grp"; c.f. sim.pwr.plot()
##' 
##' @param predName  leave as default "pvals"; c.f. sim.pwr.plot()
##' 
#' @param legend  whether to put a legend on a single plot
#' 
#' @param legend.just  the legend.justification argument in ggplot2 controling the position of the legend
#' 
#' @param x.limit  specify the plotting range for the x-axis (TPR), for better visualization
#' 
#' @param panel.label  the label on the lower-right corner of each plot (for journal figure legends)
#' 
#' @param title  the plot title for a single panel
#' 
#' @return a list of a ggplot2 object, relevant information and F-measures. 
#' 
#' @keywords internal
#' 
#' @importFrom plyr llply
#' @importFrom plyr ldply
#' 
#' @author Gu Mi
#' 
prec.plot = function(test.data.list, beta=beta, groupName = "grp", predName = "pvals", 
                     legend=TRUE, legend.just=c(0,1), x.limit=c(0,1), panel.label="A", title = "Plot Title") {
  
  plotdata = llply(test.data.list, function(x) 
    with(x, pwrdata(grp = eval(parse(text = groupName)), pred = eval(parse(text = predName)), type="PREC", beta=beta) ))
  plotdata2 = list(pwr = ldply(plotdata, function(x) x$pwr))
  names(plotdata2$pwr) = c("model.names", "recall", "precision", "pvalue", "qvalue", "Fmeasure")  
  # by default, ".id" is used for the column of model names
  
  if (legend) {
    p = 
      ggplot(plotdata2$pwr, aes(x = recall, y = precision)) +
      geom_line(size=0.5, aes(linetype=factor(model.names), colour=factor(model.names))) +
      theme_bw() +
      scale_x_continuous("", limits=x.limit) +
      scale_y_continuous("", limits=c(0,1)) +
      # in order to get legend of lines having both colour and linetype, they must have the same title and labels, etc.
      # o.w. there will be duplicate legends produced
      scale_linetype_discrete(names(test.data.list),
                              breaks = names(test.data.list),
                              labels = names(test.data.list)) +
      scale_colour_brewer(palette="Dark2", 
                          breaks = names(test.data.list),
                          labels = names(test.data.list)) +
      ggtitle(title) + 
      annotate("text", x = x.limit[2], y = 0, label=panel.label, size = 10, face="bold") +
      theme(plot.title = element_text(face="bold", size=18), 
            axis.title.x = element_text(face="bold", size=16),
            axis.title.y = element_text(face="bold", size=16, angle=90),
            axis.text.x  = element_text(size=16),
            axis.text.y  = element_text(size=16),
            plot.margin=unit(x=c(0.8, 0.8, 1, 1), units="cm"),  # top, right, bottom, and left
            legend.key.width = unit(3, "line"),
            legend.justification = legend.just,
            legend.position = legend.just,
            legend.title = element_blank(),
            legend.key = element_blank(),
            legend.text = element_text(size=16)
      )
  }
  #
  if (!legend) {
    p = 
      ggplot(plotdata2$pwr, aes(x = recall, y = precision)) +
      geom_line(size=0.5, aes(linetype=factor(model.names), colour=factor(model.names))) +
      theme_bw() +
      scale_x_continuous("", limits=x.limit) +
      scale_y_continuous("", limits=c(0,1)) +
      # in order to get legend of lines having both colour and linetype, they must have the same title and labels, etc.
      # o.w. there will be duplicate legends produced
      scale_linetype_discrete(names(test.data.list),
                              breaks = names(test.data.list),
                              labels = names(test.data.list)) +
      scale_colour_brewer(palette="Dark2", 
                          breaks = names(test.data.list),
                          labels = names(test.data.list)) +
      ggtitle(title) + 
      annotate("text", x = x.limit[2], y = 0, label=panel.label, size = 10, face="bold") +
      theme(plot.title = element_text(face="bold", size=18), 
            axis.title.x = element_text(face="bold", size=16),
            axis.title.y = element_text(face="bold", size=16, angle=90),
            axis.text.x  = element_text(size=16),
            axis.text.y  = element_text(size=16),
            plot.margin=unit(x=c(0.8, 0.8, 1, 1), units="cm"),  # top, right, bottom, and left
            legend.position = "none"
      )
  }
  return(list(pwr = p, info = plotdata2$pwr))
  
}

