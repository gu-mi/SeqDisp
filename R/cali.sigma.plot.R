

##' @title Plotting function for quantifying noise with calibration interval in ggplot2
##' 
##' @param real.data.sigma.obj  an object from \code{quantify.noise} where the second evaluate is TRUE
##' @param sim.sigma.obj  an object from \code{quantify.noise} where the third evaluate is TRUE
##' @param sigma.v  a vector of true sigma (to be plotted on the x-axis)
##' @param cali.sigma  the calibrated sigma estimate (need to see the plot and then determine this value)
##' @param left.ci  the left bound for the calibration interval (need to see the plot and then determine this value)
##' @param right.ci  the right bound for the calibration interval (need to see the plot and then determine this value)
##' @param title  main title of the plot
##' 
##' @return a ggplot2 object
##' 
##' @export
##' 
##' @author Gu Mi
##' 
cali.sigma.plot = function(real.data.sigma.obj, 
                           sim.sigma.obj, 
                           sigma.v = c(0.5, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.5),
                           cali.sigma = 1.0,
                           left.ci = 0.8,
                           right.ci = 1.2,
                           title = NULL){  
  
  sigma = real.data.sigma.obj$sigma.v  # true sigma values
  sigmasq = sigma^2
  sigmas = unique(sim.sigma.obj$sim.sigma.v.medians)
  myFit = lm(sigmas ~ sigma + sigmasq)
  myPredict = as.data.frame(predict(myFit, interval="prediction"))
  myData = cbind(myPredict, sigma)
  
  # begin plotting
  
  sigmap = ggplot(data=sim.sigma.obj, aes(x = sigma.v, y = sim.sigma.v.medians)) + 
    
    # all sigma's:
    geom_point(data=sim.sigma.obj, aes(x = sigma, y = sigma.hat), colour="Gray", size=1, shape=1) + 
    # median of each of the unique-level sigma's:
    geom_point(data=sim.sigma.obj, aes(x = sigma.v, y = sim.sigma.v.medians), colour="Blue", size=1.5, shape=20) + 
    # horizontal line:
    geom_line(data=real.data.sigma.obj, aes(x = sigma.v, y = real.data.sigmas), linetype="solid", colour="Black", size=0.3) + 
    # fitted line:
    geom_line(data=myData, x = sigma, y = myPredict[,1], colour="Black", linetype="solid", size=0.3) + 
    # lwr:
    geom_line(data=myData, x = sigma, y = myPredict[,2], colour="Black", linetype="dashed", size=0.2) + 
    # upr:
    geom_line(data=myData, x = sigma, y = myPredict[,3], colour="Black", linetype="dashed", size=0.2) + 
    # THE point:
    geom_point(data=real.data.sigma.obj, x = cali.sigma, aes(y = real.data.sigmas[1]), size = 0.6, colour = "Red", shape=21) +
    # calibrations:
    geom_vline(xintercept = left.ci, colour="Black", linetype = "longdash", size=0.2) +
    geom_vline(xintercept = right.ci, colour="Black", linetype = "longdash", size=0.2) +
    #
    scale_x_continuous(expression(bold(paste("True ", sigma, sep=""))), breaks = sigma.v) +
    scale_y_continuous(expression(bold(paste("Estimated ", hat(sigma), sep=""))), breaks = sigma.v, 
                       limits = c(min(sim.sigma.obj$sim.sigma.v.medians), max(sim.sigma.obj$sim.sigma.v.medians))) +
    ggtitle(title) + 
    theme_bw() + 
    theme(plot.title = element_text(face="bold", size=8),
          # axis labels
          axis.title.x = element_text(face="bold", size=8),
          axis.title.y = element_text(face="bold", size=8, angle=90),
          # axis tick labels
          axis.text.x = element_text(angle=0, vjust=1, size=6), 
          axis.text.y = element_text(size=6),
          legend.position="none"
    )
  return(sigmap)
}