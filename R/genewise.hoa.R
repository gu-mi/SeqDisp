

##' @title Implement genewise dispersion approach with LRT-HOA test
##' 
##' @param nb.data
##' @param grp.ids
##' @param grp1
##' @param grp2
##' @param R
##' @param print.level
##' 
##' @return a data frame where each row is from hoa.1u() function output
##' 
##' @export
##' 
##' @author Yanming Di
##' 
genewise.hoa = function(nb.data, grp.ids, grp1, grp2, R = 100, print.level=1) {
  
  counts = nb.data$counts;
  eff.lib.sizes = nb.data$eff.lib.sizes;
  
  ## Extract relevant columns from the counts matrix
  id1 = grp.ids %in% grp1;
  id2 = grp.ids %in% grp2;
  
  if (print.level>0) {
    message(sprintf("DE test betwee group %s and group %s.\n", grp1, grp2));
  }
  
  y = cbind(counts[, id1, drop=FALSE], counts[, id2, drop=FALSE]);
  
  s = c(eff.lib.sizes[id1], eff.lib.sizes[id2]);
  
  ## Construct the model matrix
  treatment = unlist(list(grp.ids[id1], grp.ids[id2]));
  x = model.matrix(~factor(treatment, levels = c(grp1, grp2)));
  colnames(x) = c(grp1, grp2);
  
  ## The null hypothesis
  beta0 = c(NA, 0);
  
  ## Test DE using HOA
  ## set.seed(999);
  
  m = nrow(y);
  n = ncol(y);
  ## m = 100;
  
  res = data.frame(
    ## statistic = I(matrix(NA, m, 3)),
    ## p.value = I(matrix(NA, m, 3)),
    t.hoa=numeric(m), p.hoa=numeric(m),
    t.lr=numeric(m), p.lr=numeric(m),
    t.wald=numeric(m), p.wald=numeric(m),
    kappa.hat = numeric(m), beta.hat = I(matrix(NA, m,2)), mu.hat = I(matrix(NA,m,n)),
    kappa.tilde = numeric(m), beta.tilde = I(matrix(NA, m, 2)), mu.tilde = I(matrix(NA, m,n)),
    zsim = numeric(m),
    alternative=character(m),
    msg=character(m),
    stringsAsFactors=FALSE);
  
  ## system.time({
  for (i in 1:m) {
    res[i,] = hoa.1u(y[i,], s = s, x = x, beta0 = beta0, R=R);
  }
  
  res;
}