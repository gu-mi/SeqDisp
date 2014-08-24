

##' @title Fit polynomial gamma log-linear regression with linear, quadratic or cubic terms
##' 
##' @param counts  read counts matrix
##' @param x  design matrix
##' @param model  choose from "linear", "quadratic" or "cubic"
##' 
##' @return a list showing the significance, R^2 of each of the fitted polynomial regression
##' 
##' @export
##' 
##' @author Gu Mi
##' 
poly.gamma.reg = function(counts, x, model = "linear"){
  
  m = dim(counts)[1]
  n = dim(counts)[2]
  
  stopifnot(model %in% c("linear", "quadratic", "cubic"))
  
  mu.mom = expandAsMatrix(rowMeans(counts), dim=c(m,n))
  re.freq = (mu.mom / (matrix(1, m, 1) %*% matrix(colSums(counts), 1, n)))[ ,1] 
  #qt.re.freq = quantile(re.freq, c(0.001, 0.999))
  phi.mom = (rowSums((counts - mu.mom)^2) - rowSums(mu.mom))/rowSums(mu.mom^2)  # may have NaN
  #id = (phi.mom > 0 & !is.nan(phi.mom) & re.freq > qt.re.freq[1] & re.freq < qt.re.freq[2])  
  id = (phi.mom > 0 & !is.nan(phi.mom))
  # may discard some phi.hat here not in the plotting: this "id" is used "globally" to subset
  n.pt = sum(id)
  
  phi.hat = phi.mom
  
  if (model == "linear"){
    # linear trend (line)
    s1 = glm(phi.hat[id] ~ poly(log(re.freq[id]), degree=1), family=Gamma(link="log"))
    pred.resp.1 = predict(s1, 
                          type = "response", 
                          newdata = data.frame(cbind(rep(1,length(re.freq[id])), log(re.freq[id]))))
    
    s0 = glm(phi.hat[id] ~ 1, family=Gamma(link="log"))
    SS1 = sum((log(phi.hat[id]) - log(fitted.values(s0)))^2)   
    SS2 = sum((log(phi.hat[id]) - log((pred.resp.1)))^2)    
    R.sq = (SS1 - SS2)/SS1  
    
    res = list(sig = summary(s1)$coef[,"Pr(>|t|)"],
               R.sq = R.sq)
  }
  
  if (model == "quadratic"){
    # quadratic trend (curve)
    s2 = glm(phi.hat[id] ~ poly(log(re.freq[id]), degree=2), family=Gamma(link="log"))
    pred.resp.2 = predict(s2, 
                          type = "response", 
                          newdata = data.frame(cbind(rep(1,length(re.freq[id])), log(re.freq[id]), (log(re.freq[id]))^2)))
    
    s0 = glm(phi.hat[id] ~ 1, family=Gamma(link="log"))
    SS1 = sum((log(phi.hat[id]) - log(fitted.values(s0)))^2)   
    SS2 = sum((log(phi.hat[id]) - log((pred.resp.2)))^2)    
    R.sq = (SS1 - SS2)/SS1  
    
    res = list(sig = summary(s2)$coef[,"Pr(>|t|)"],
               R.sq = R.sq)
  }
  
  if (model == "cubic"){
    # cubic trend (curve)
    s3 = glm(phi.hat[id] ~ poly(log(re.freq[id]), degree=3), family=Gamma(link="log"))
    pred.resp.3 = predict(s3, 
                          type = "response", 
                          newdata = data.frame(cbind(rep(1,length(re.freq[id])), log(re.freq[id]), 
                                                     (log(re.freq[id]))^2, (log(re.freq[id]))^3)))      
    
    s0 = glm(phi.hat[id] ~ 1, family=Gamma(link="log"))
    SS1 = sum((log(phi.hat[id]) - log(fitted.values(s0)))^2)   
    SS2 = sum((log(phi.hat[id]) - log((pred.resp.3)))^2)    
    R.sq = (SS1 - SS2)/SS1  
    
    res = list(sig = summary(s3)$coef[,"Pr(>|t|)"],
               R.sq = R.sq)
  }
  
  return(res)
}

