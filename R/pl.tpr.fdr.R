
#' @title  Plot True Positive Rate (TPR) vs. False Discovery Rate (FDR) using ggplot2
#' 
#' @details This function calls \code{sim.data.prepare} to generate simulated read counts with different dispersion models,
#' varying percentages of DE genes, log fold changes and random normal noise levels. This function calls \code{sim.pwr.plot}
#' to generate a ggplot2 object for a single scenario (DE percentage, log FC and noise level combination) and then loop over
#' all considered scearios (panels of plots). This function also output .txt files (for each panel) the TPR, FDR, p-value and
#' q-values for each of the models considered (as specified in the \code{model.vec} argument).
#' 
#' @param dt  the full dataset to be loaded and used as a template for the simulation
#' 
#' @param sub.col  column indices used to subset the full dataset
#' 
#' @param grp.ids  group ids to distinguish between the two groups
#' 
#' @param disp.vec  the dispersion types to be simulated. Specify "Linear", "Quadratic" or "Non-parametric"
#' 
#' @param init.fit  (logical) whether to estimate the dispersion model and save as a template (for later use)
#' 
#' @param m  number of genes (default = 5000 genes) subsetted for actual analysis
#' 
#' @param n  number of samples (default = 6 samples)
#' 
#' @param s  library size (default = 1e7, 10 million reads per library)
#' 
#' @param de.vec  a vector of percentages of DE genes (default = 0.2)
#' 
#' @param fc.vec  a vector of log fold-changes, i.e. the value of beta2 (default = log(3.0), for FC = 3.0).
#' If a vector of fold changes, will call \code{sim.data.prepare}; 
#' If NULL, will call \code{sim.data.prepare2}.
#' 
#' @param noise.vec  a vector of random normal(0, sigma) noise to phi: phi * exp(rnorm(m, 0, sigma)). 
#' If any of the noise.vec elements equals 0 (default), data without noise are simulated
#' 
#' @param legend.just  the legend.justification argument in ggplot2 controling the position of the legend.
#' 
#' @param x.limit  specify the plotting range for the x-axis (TPR), for better visualization
#' 
#' @param seed  random generator seed for reproducible research
#' 
#' @param out.path  directory path where the intermediate outputs are stored (and then fetched)
#' 
#' @param model.vec  a vector of NB dispersion models to be compared for power-robustness, including: 
#' "Common", "NBP", "NBQ", "NBS", "Trended", "Genewise", "Tagwise-Common", "Tagwise-Trend", "QL", "QLShrink" and "QLSpline"
#' 
#' @return a list of two component: (1) information for ggplot2 plotting; (2) useful quantities in calculating TPR and FDR.
#' Also, this function output .txt files for future use.
#' 
#' @export
#'
#' @author Gu Mi 
#' 
pl.tpr.fdr = function(dt = NULL, 
                      sub.col = 1:6,
                      grp.ids = c(rep(1,length(sub.col)/2), rep(2,length(sub.col)/2)),
                      disp.vec = NULL,
                      init.fit = TRUE,
                      m = 5000, 
                      n = 6,
                      s = 1e7, 
                      de.vec = 0.2,
                      fc.vec = log(3.0),
                      noise.vec = 0,
                      legend.just = c(0,1),
                      x.limit = c(0,1),
                      seed = 539, 
                      out.path = getwd(),
                      model.vec = NULL) {
  
  n.rows = length(noise.vec)
  m.pwr = rep( list(NA), n.rows) 
  legend = c(TRUE, rep(FALSE, (n.rows-1)))  # put legend only on the 1st panel
  panel.label = LETTERS[1:n.rows]
  
  for ( i in seq_len(length(noise.vec)) ){
    
    if (!is.null(fc.vec)) {
      # sim.data.prepare(): use "MAPL" for the dispersion estimation
      d.prep = sim.data.prepare(dt = dt, sub.col = sub.col, grp.ids = grp.ids, init.fit = init.fit, 
                                disp = disp.vec[i], m = m, n = n, s = s, perc.de = de.vec[i], log.fc = log(fc.vec[i]), sigma = noise.vec[i], 
                                seed = seed, out.path = out.path)
      # sim.pwr.plot()      
      p.info = sim.pwr.plot(sim.data.obj = d.prep,
                            type = "TPR-FDR",
                            model.vec = model.vec, 
                            legend=legend[i],
                            legend.just = legend.just,
                            x.limit = x.limit,
                            panel.label=panel.label[i],
                            title.note = bquote(
                              paste(.(disp.vec[i]), " Trend: ", sigma, " = ", .(noise.vec[i]), 
                                    "; DE% = ", .(de.vec[i]*100), "%", sep=""))
                            #title.note = bquote(paste(.(disp.vec[i]), " Trend with ", .(de.vec[i]*100), "% DE genes", sep=""))
      )
      m.pwr[[i]] = p.info$pwr
      info.table = p.info$info
      colnames(info.table) = c("model", "TPR", "FDR", "pvalue", "qvalue")
      write.table(x=info.table, quote=FALSE, row.names=FALSE, col.names=TRUE, 
                  file=file.path(out.path, paste0("TF-panel-", panel.label[i] , "-DE", de.vec[i], "-FC", fc.vec[i], ".txt")))
    }
    #
    if (is.null(fc.vec)) {
      # sim.data.prepare2(): use "MAPL" for the dispersion estimation
      d.prep = sim.data.prepare2(dt = dt, sub.col = sub.col, grp.ids = grp.ids, init.fit = init.fit, 
                                 disp = disp.vec[i], m = m, n = n, s = s, perc.de = de.vec[i], sigma = noise.vec[i], 
                                 seed = seed, out.path = out.path)
      # sim.pwr.plot()      
      p.info = sim.pwr.plot(sim.data.obj = d.prep,
                            type = "TPR-FDR",
                            model.vec = model.vec, 
                            legend=legend[i],
                            legend.just = legend.just,
                            x.limit = x.limit,
                            panel.label=panel.label[i],
                            title.note = bquote(
                              paste(.(disp.vec[i]), " Trend: ", sigma, " = ", .(noise.vec[i]), 
                                    "; DE% = ", .(de.vec[i]*100), "%", sep=""))
                            #title.note = bquote(paste(.(disp.vec[i]), " Trend with ", .(de.vec[i]*100), "% DE genes", sep=""))
      )
      m.pwr[[i]] = p.info$pwr
      info.table = p.info$info
      colnames(info.table) = c("model", "TPR", "FDR", "pvalue", "qvalue")
      write.table(x=info.table, quote=FALSE, row.names=FALSE, col.names=TRUE, 
                  file=file.path(out.path, paste0("TF-panel-", panel.label[i] , "-DE", de.vec[i], "-Sigma", noise.vec[i], ".txt")))
    }
    
  }
  
  return(list(m.pwr=m.pwr))
  
}