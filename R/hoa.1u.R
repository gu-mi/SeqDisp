
##' HOA test for regresssion coefficeint in a single NB regression model
##'
##' @title HOA test for regresssion coefficeint in a NB regression model
##' @param y 
##' @param s 
##' @param x 
##' @param beta0 
##' @param alternative 
##' @param R 
##' @param print.level 
##' @return a data frame
##' @author Yanming Di
##' 
##' @keywords internal
##' 
hoa.1u = function(y, s, x, beta0,
                  alternative = "two.sided",
                  R=100,
                  print.level=0) {
  
  if (print.level > 0) 
    message("HOA test for regression coefficients in a NB regression model with unknown dispersion.");
  
  ## Initialize test statistic and p-values.
  ## If the funtion terminated abnormally, these values will be returned.
  r = rsim = w = NA;
  p = psim = pstar = p.wald =  NA;
  kappa.hat = NA;
  beta.hat = NA;
  mu.hat = NA;
  kappa.tilde = NA;
  beta.tilde = NA;
  mu.tilde = NA;
  zsim = NA;
  msg = "";
  
  alternatives = c("two.sided", "less", "greater");
  alt = pmatch(alternative, alternatives);
  if (is.na(alt)) {
    stop("alternative should be \"two.sided\", \"less\" or \"greater\"")
  }
  alternative = alternatives[alt];
  
  ## Prepare return values
  return.value = function(){
    if (alt == 1) {
      p = min(p * 2, 1)
      p.wald = min(p.wald * 2, 1)
      psim = min(psim * 2, 1);
      ## p.score = min(p.score * 2, 1)
    }
    
    ## statistic = matrix(c(rsim, r, w), 1, 3);
    ## p.value = matrix(c(psim, p, p.wald), 1, 3);
    ## names(statistic) = c("HOA", "LR", "Wald");
    ## names(p.value) = c("HOA", "LR", "Wald");
    
    data.frame(
      ## statistic = I(statistic),
      ## p.value = I(p.value),
      t.hoa = rsim, p.hoa = psim,
      t.lr = r, p.lr = p,
      t.wald = w, p.wald = p.wald,
      ## r=r, p=p, rsim=rsim, psim=psim, w=w, p.wald=p.wald,
      kappa.hat = kappa.hat,
      beta.hat = I(matrix(beta.hat, 1, length(beta.hat))),
      mu.hat = I(matrix(mu.hat, 1, length(mu.hat))),
      kappa.tilde = kappa.tilde,
      beta.tilde = I(matrix(beta.tilde, 1, length(beta.tilde))),
      mu.tilde = I(matrix(mu.tilde, 1, length(mu.tilde))),
      zsim = zsim,
      alternative=alternative,
      msg=msg,
      stringsAsFactors=FALSE);
  };
  
  if (any(y<0)) stop("some y are negative!")
  if (all(y<=0)) return(return.value());
  
  ## Indices to nuisance parameters and to the parameter to be tested
  idh = !is.na(beta0);
  nh = sum(idh);
  if (sum(idh) != 1) 
    stop("The hypothesis should be one dimensional.")
  idn = is.na(beta0)
  nn = sum(idn)
  
  idn1 = c(TRUE, idn);
  idh1 = c(FALSE, idh);
  
  ## beta0 = c(NA, 0);
  n = length(y);
  n.pars = length(beta0) + 1;
  
  ## Fit full model
  full = fit.nb.regression.1(y, s, x);
  kappa.hat = full$kappa;
  beta.hat = full$beta;
  mu.hat = full$mu;
  l.hat = full$l;
  j.hat = full$j;
  
  ## Fit redueced model
  reduced = fit.nb.regression.1(y, s, x, beta0);
  kappa.tilde = reduced$kappa;
  beta.tilde = reduced$beta;
  mu.tilde = reduced$mu;
  l.tilde = reduced$l;
  j.tilde = reduced$j;
  
  ## Determine 
  if (alt == 2 | (alt == 1 & beta.hat[idh] < beta0[idh])) {
    lower.tail = TRUE
  } else {
    lower.tail = FALSE
  }
  
  ## Special care is needed if kappa.hat is large
  ## Reduce to Poisson regression
  if (kappa.hat > 1e8) {
    msg = "[kappa.hat>1e8]";
    big.kappa.hat = TRUE;
  } else {
    big.kappa.hat = FALSE;
  }
  
  ## Wald test
  if (big.kappa.hat) {
    ## Wald test (using observed infomation matrix)
    res = try({
      (beta.hat[idh] - beta0[idh])/sqrt(solve(j.hat[-1,-1])[idh, idh]);
    });
    
    if ("try-error" %in% class(res)) {
      msg=paste(msg, "[error computing wald test]");
      w = NA;
      p.wald = NA;
    } else {
      w = res;
      p.wald = pnorm(w, lower.tail = lower.tail)
    }
  } else {
    ## Wald test (using observed infomation matrix)
    res = try({
      (beta.hat[idh] - beta0[idh])/sqrt(solve(j.hat)[idh1, idh1]);
    });
    
    if ("try-error" %in% class(res)) {
      msg=paste(msg, "[error computing wald test]");
      w = NA;
      p.wald = NA;
    } else {
      w = res;
      p.wald = pnorm(w, lower.tail = lower.tail)
    }
  }
  
  ## If l.hat and l.tilde are the same 
  if (l.hat <= l.tilde + sqrt(.Machine$double.eps)) {
    msg=paste(msg, "[l.hat<=l.tilde]");
    
    r = 0;
    rsim = 0;
    p = 0.5;
    psim = 0.5;
    
    return(return.value());
  }
  
  ## LR test
  r = sign(beta.hat[idh] - beta.tilde[idh]) * sqrt(2 * (l.hat - l.tilde));
  p = pnorm(r, lower.tail = lower.tail);
  
  ## Sample space derivatives
  S.hat = matrix(0, n.pars, n.pars);
  ## S.hat[1,1] = n /2 / sigma.sq.tilde^2;
  ## S.hat[-1,-1] =  t(x) %*% x / sigma.sq.tilde;
  ## S.hat[-1,1] = matrix((mu.hat - mu.tilde)/sigma.sq.tilde^2, 1, n) %*% x;
  
  q.hat = numeric(n.pars);
  ## q.hat[1] = - n/2/sigma.sq.hat + n/2/sigma.sq.tilde;
  ## q.hat[-1] = matrix((mu.hat - mu.tilde)/sigma.sq.tilde, 1, n) %*% x;
  
  ## z = r * det(i.hat)/sqrt(det(j.hat))/det(S.hat) * sqrt(det(j.tilde[idn1, idn1]))/(solve(S.hat) %*% q.hat)[idh1];
  ## rstar = r - (1/r) * log(z)
  ## pstar = pnorm(rstar, lower.tail = lower.tail);
  
  ll = function(par, y) {
    ## kappa = exp(par[1]);
    kappa = par[1];
    beta = par[-1];
    mu = s * exp(x %*% beta);
    ## log.likelihood.nb(kappa, mu, y);
    sum(dnbinom(y, size=kappa, mu=mu, log=TRUE));
  }
  
  D1 = function(par, y) {
    
    kappa = par[1];
    beta = par[-1];
    mu = s * exp(x %*% beta);
    v = mu + mu^2/kappa;
    
    ## log.likelihood.nb(kappa, mu, y);
    ## sum(dnbinom(y, size=kappa, mu=mu, log=TRUE));
    dl.dbeta = matrix((y - mu) * mu/v, 1, n) %*% x;
    dl.dkappa = sum(digamma(kappa + y) - digamma(kappa) + log(kappa) + 1 - (log(mu+kappa) + (y+kappa)/(mu+kappa)));
    c(dl.dkappa, dl.dbeta);
  }
  
  theta.hat = c(kappa.hat, beta.hat);
  theta.tilde = c(kappa.tilde, beta.tilde);
  
  ## R = 10;
  
  ## system.time({
  ## Simulate normal data
  ## debug(D1);
  meanAll = rep(0,2*n.pars+1);
  prodAll = matrix(0,2*n.pars+1, 2*n.pars+1);
  for (i in 1:R) {
    y0 = rnbinom(n, size=kappa.hat, mu=mu.hat);
    l1 = ll(theta.hat, y=y0); 
    l0 = ll(theta.tilde, y=y0); 
    ## uhat = grad(ll, theta.hat, y=y0);
    ## uhyp = grad(ll, theta.tilde, y=y0);
    uhat = D1(theta.hat, y=y0);
    uhyp = D1(theta.tilde, y=y0);
    uhh = c(uhat, uhyp, l1-l0);
    meanAll = meanAll + uhh / R;
    prodAll = prodAll +  tcrossprod(uhh) / R;
  }
  
  covAll = prodAll * R/(R-1) - tcrossprod(meanAll)  * R/(R-1);
  
  S.sim =  covAll[1:n.pars, (n.pars+1):(2*n.pars)];
  i.sim = covAll[1:n.pars, 1:n.pars];
  ## iHatInv <- solve(iHat)
  q.sim  =  covAll[1:n.pars, 2*n.pars+1];
  
  ## HOA test
  if (big.kappa.hat) {
    ## Recude to Poisson regression
    ## Remove the kappa components from the likelihood quantities
    i.sim = i.sim[-1,-1];
    j.hat = j.hat[-1,-1];
    S.sim = S.sim[-1,-1];
    j.tilde = j.tilde[-1,-1];
    q.sim = q.sim[-1];
    
    res = try({
      r * det(i.sim)/sqrt(det(j.hat))/det(S.sim) * sqrt(det(j.tilde[idn, idn, drop=FALSE]))/(solve(S.sim) %*% q.sim)[idh];
    });
    
    if ("try-error" %in% class(res) | is.na(res)) {
      msg=paste(msg, "--error computing zsim--");
      
      message("kappa.hat=", kappa.hat,"\n");
      zsim = NA;
      rsim = NA;
      psim = NA;
    } else {
      if (res<0) {
        msg = paste(msg, "zsim<0");
        zsim = 1;
      } else {
        zsim = res;
      }
      rsim = r - (1/r) * log(zsim);
      psim = pnorm(rsim, lower.tail = lower.tail);
    }
    
  } else {
    res = try({
      r * det(i.sim)/sqrt(det(j.hat))/det(S.sim) * sqrt(det(j.tilde[idn1, idn1, drop=FALSE]))/(solve(S.sim) %*% q.sim)[idh1];
    });
    
    if ("try-error" %in% class(res) | is.na(res)) {
      msg=paste(msg, "[error computing zsim]");
      
      message("kappa.hat=", kappa.hat,"\n");
      zsim = NA;
      rsim = NA;
      psim = NA;
    } else {
      
      if (res<0) {
        msg = paste(msg, "[zsim<0]");
        zsim = 1;
      } else {
        zsim = res;
      }
      rsim = r - (1/r) * log(zsim);
      psim = pnorm(rsim, lower.tail = lower.tail);
    }
    
  }
  ## });
  
  return.value();
}
