#include <RcppArmadillo.h>
using namespace Rcpp;

//[[Rcpp::depends(RcppArmadillo)]]


void softthresh(arma::vec x, double thresh){
  arma::vec::iterator it = x.begin();
  arma::vec::iterator it_end = x.end();
  int sgn;
  double shrink;
  for(;it != it_end; ++it){
    sgn = arma::sign(*it);
    shrink = (*it)*sgn - thresh;
    shrink = (shrink >= 0) ? shrink : 0;
    *it = sgn*shrink;
  }
}

// [[Rcpp::export]]
List ls_f_sparse(int bn, arma::uvec bs, arma::uvec ix,arma::uvec iy,arma::vec gam,
                 int nobs, int nvars, arma::mat const& x, arma::vec const &y,
                 arma::vec pf, int dfmax, int pmax, int nlam, double flmin,
                 arma::vec ulam, double eps, int maxit, int nalam,
                 int npass, int jerr,
                 double alsparse){
  // double big=9.9E30;
  double mfl = 1.0E-6;
  int mnlam = 6;
  int mnl;
  
  // - - - local declarations - - -
  double max_gam, dif, al, alf;// sg, d;
  arma::vec b(nvars);
  arma::vec oldbeta(nvars);
  arma::vec r = y; // Residue y-beta_k*x etc
  arma::vec oldb;
  arma::vec dd;
  arma::uvec idx(pmax); 
  arma::vec oidx(bn);
  int g, j, l, ni, me, startix, endix;
  // - - - Aaron's declarations
  double snorm;
  arma::vec t_for_s(bn); // this is for now just 1/gamma
  double tea; // this takes the place of 't' in the update step for ls
  arma::vec s; // takes the place of 'u' in update for ls
  // int vl_iter; // for iterating over columns(?) of x*r
  // int kill_count=0;
               
  double tlam, lama, lam1ma, al0=0.0, testr=0.0, tmp;
  int jx;
  arma::uvec jxx(bn);
  arma::vec ga(bn);
  max_gam = gam.max();
  
  if(max(pf) <= 0.0){
    jerr=10000; // Fatal error, we quit
    List out = List::create(Named("jerr")=jerr);
    return out;
  }
  pf = arma::clamp(pf, 0.0, pf.max());
  jxx.zeros();
  al = 0.0;
  mnl = (nlam >= mnlam) ? mnlam : nlam; // min
  // ----- output in cpp -----
  // arma::vec b0(nlam);
  arma::sp_mat beta(nvars, nlam);
  arma::vec alam(nlam);
  arma::uvec nbeta(nlam);
  
  b.zeros();
  oldbeta.zeros();
  idx.zeros();
  oidx.zeros();
  npass = 0;
  ni = 0;
  alf = 0.0;
  t_for_s = 1/gam;
  // --------- lambda loop ----------------------------
  if(flmin < 1.0){ // THIS is the default...
    flmin = (mfl >= flmin) ? mfl : flmin; // just sets a threshold above zero 
    alf = pow(flmin, 1.0/(nlam-1.0));
  }
  // Rcout << alf << std::endl;
  arma::vec vl = x.t() * r / nobs; // Note r gets updated in middle and inner loop 
  al0 = 0.0;
  // For each group...      
  for(g = 0; g<bn; g++) ga(g) = arma::norm(vl.subvec(ix(g),iy(g)));
  al0 = arma::norm(vl, "inf"); // Infty norm of X'y, big overkill for lam_max
  // Rcout << al0 << std::endl;
  al = al0; // ! this value ensures all betas are 0
  l = -1;
  tlam = 0.0;
  while (l < nlam-1){ // This is the start of the loop over all lambda values...
    al0 = al; // store old al value on subsequent loops, first set to al
    if (flmin>=1.0){ // user supplied lambda value, break out of everything
      l++;
      al=ulam(l);
    } else {
      if (l > 0){ // have some active groups
	      al*=alf;
        tmp = 2.0*al-al0;
	      tlam = (tmp >= 0) ? tmp : 0.0; // Here is the strong rule...
	      l++;
      } else if (l < 0){ //Trying to find an active group
	      al*=.99;
	      tlam=al;
      }
    }
    // Rcout << "l = " << l << std::endl;
    lama = al*alsparse;
    lam1ma = al*(1-alsparse);
    // This is the start of the algorithm, for a given lambda...
    for(g = 0; g < bn; g++){ 
      if(jxx(g) == 1)  continue;
      if(ga(g) > pf(g)*tlam*(1-alsparse)) jxx(g)++; // Implementing the strong rule
    }
    // ! --------- outer loop ---------------------------- ! 
    do { // while condition checked way down below (jx==1)
      // Rcout << "This is the outer loop." << std::endl;
      if(ni>0){
	      for(j=0; j<ni; j++){
	        g=idx(j);
	        oldbeta.subvec(ix(g),iy(g))=b.subvec(ix(g),iy(g));
	      }
      }
      // --middle loop-------------------------------------
      for(;;) { // loops forever, we rely on one of the break statements
        // Rcout << "This is the middle loop." << std::endl;
	      npass++;
	      dif=0.0;
	      for(g=0; g<bn; g++){
	        if(jxx(g) == 0) continue;
	        startix=ix(g);
	        endix=iy(g);
	        oldb=b.subvec(startix, endix);
	        s = x.cols(startix,endix).t() * r / nobs;
	        s*=t_for_s(g);
	        s += oldb;
	        softthresh(s, lama*t_for_s(g));
	        snorm = arma::norm(s);
	        tea = snorm - t_for_s(g)*lam1ma*pf(g);
	        if (tea>0.0){
	          b.subvec(startix,endix) = s*tea/snorm;
	        } else {
	          b.subvec(startix,endix).zeros();
	        } 
	        dd=b.subvec(startix, endix)-oldb;
	        if(arma::any(dd)){ //seeming issue with .is_zero()
	          tmp = gam(g)*gam(g)*arma::dot(dd,dd);
	          dif = (dif >= tmp) ? dif : tmp;
	          r -= x.cols(startix, endix) * dd;
	          if(oidx(g)==0){// Here is where middle loop is different;
	            // if group g was not in oidx (active), and the
	            // difference was nonzero, put it in active (ni)
	            ni++;
	            if(ni>pmax) break;
	            oidx(g)=ni;
	            idx(ni-1)=g;
	          }
	        }
	      }
	      if (ni > pmax) break;
	      if (dif < eps) break;
	      if (npass > maxit){ // too many iterations, but we return everything
	        jerr = -l+1;
	        beta.resize(nvars,l-1);
	        alam.resize(l-1);
	        List out = List::create(
	          Named("jerr") = jerr,
            Named("nbeta") = nbeta,
            Named("beta") = beta,
            Named("alam") = alam,
            Named("nalam") = nalam,
            Named("npasses") = npass,
            Named("nalam") = l-1);
	        return out;
	      }
	      // --inner loop----------------------
	      do { // do...while() loop, see the end for the condition
	        // Rcout << "This is the inner loop." << std::endl;
	        npass++;
	        dif=0.0;
	        for(j=0; j<ni; j++){
	          g=idx(j);
	          startix=ix(g);
	          endix=iy(g);
	          oldb=b.subvec(startix, endix);
	          s = x.cols(startix,endix).t() * r / nobs;
	          s *= t_for_s(g);
	          s += oldb;
	          softthresh(s, lama*t_for_s(g));
	          snorm = arma::norm(s);
	          tea = snorm - t_for_s(g)*lam1ma*pf(g);
	          if (tea>0.0){
	            b.subvec(startix,endix) = s*tea/snorm;
	          } else {
	            b.subvec(startix,endix).zeros();
	          }
	          dd=b.subvec(startix, endix)-oldb;
	          if(!dd.is_zero()){
	            tmp = gam(g)*gam(g)*arma::dot(dd,dd);
	            dif = (dif >= tmp) ? dif : tmp;
	            r -= x.cols(startix, endix) * dd;
	          }
	        }
	        if (npass > maxit){ // too many iterations, but we return everything
	          jerr = -l+1;
	          beta.resize(nvars,l-1);
	          alam.resize(l-1);
	          List out = List::create(
	            Named("jerr") = jerr,
	            Named("nbeta") = nbeta,
	            Named("beta") = beta,
	            Named("alam") = alam,
	            Named("npasses") = npass,
	            Named("nalam") = alam.n_elem);
	          return out;
	        }
	      } while(dif > eps);  // End inner loop
      } // End middle loop                    
      if(ni>pmax) break;
      // !--- final check ---------- ! This checks which violate KKT condition
      // Rcout << "This is the KKT check." << std::endl;
      jx = 0;
      for(j=0; j<nvars; j++){
        testr = max_gam*(b(j)-oldbeta(j))/(1+std::fabs(b(j)));
        testr *= testr;
        if(testr >= eps){
          jx = 1;
          break;
        }
      }
      if (jx == 1){ 
        // Rcout << "testr = " << testr << std::endl;
        continue; // return to top of outer loop
      }
      vl = x.t() * r / nobs;
      for(g=0; g<bn; g++){
        if(jxx(g) == 1) continue;
        startix=ix(g);
        endix=iy(g);
        s = x.cols(startix,endix).t() * r / nobs;
        softthresh(s, lama);
        snorm = arma::norm(s);
        ga(g) = snorm;
        if (ga(g) > pf(g)*lam1ma){
          jxx(g) = 1;
          jx = 1;
        }
      }
      if (jx == 1){ 
        // Rcout << jxx << std::endl;
        continue; // return to top of outer loop
      }
    } while(jx==1); // Ends outer loop
    // ---------- final update variable and save results------------
    if (l < 0){
      if(!arma::any(jxx)){
        continue; // don't save anything, we're still decrementing lambda
      } else {
        l=1;
        tmp = (alf >= .99) ? alf : .99;
        alam(0) = al / tmp; // store previous, larger value
      }
    }
    // Rcout << "This is us saving output." << std::endl;
    // Rcout << "ni = " << ni << std::endl;
    if (ni > 0) { // save things first, not sure why theirs throws it out
      me=0;
      for(j=0; j<ni; j++){
        g=idx(j);
        if(arma::any(b.subvec(ix(g),iy(g)))){
          beta(arma::span(ix(g),iy(g)), l)=b.subvec(ix(g),iy(g));
          // Rcout << beta(arma::span(ix(g),iy(g)), l) << std::endl;
          me++;
        }
      }
    }
    nbeta(l)=ni;
    alam(l)=al;
    if (l < mnl) continue;
    if (ni>pmax || me>dfmax) { // now check if too many vars or groups
      jerr = -10000-l;
      break;
    }
  } // end lambda loop
  if(l < nlam-1){ // stopped early for maxit or 
    beta.resize(nvars,l);
    alam.resize(l);
  }
  List out = List::create(
    Named("jerr") = jerr,
    Named("nbeta") = nbeta,
    Named("beta") = beta,
    Named("alam") = alam,
    Named("nalam") = alam.n_elem,
    Named("npasses") = npass);
  return out;
}


/*** R
ls_sparse_cpp <- function(bn, bs, ix, iy, nobs, nvars, x, y, pf, dfmax, 
                      pmax, nlam, flmin, ulam, eps, maxit, vnames, group, intr, 
                      asparse, standardize) {
  #################################################################################
  # call Fortran core
  intercept = 0L
  if(intr){
    ym = mean(y)
    xm = colMeans(x)
    x = sweep(x,2,xm)
    y = y-ym
  }
  if(standardize){
    xs = sqrt(colSums(x^2))
    x = sweep(x,2,xs,"/")
  }
  gamma <- rep(NA, bn)
  for (g in 1:bn) gamma[g] <- RSpectra::svds(x[,ix[g]:iy[g]],1,0,0)$d^2
  gamma <- gamma/nobs
  gamma <- as.double(gamma)
  fit <- ls_f_sparse(bn, bs, ix, iy, gamma, nobs, nvars, as.matrix(x), 
                  as.matrix(y, ncol=1), pf, dfmax, pmax, nlam, flmin, 
                  ulam, eps, maxit, #intercept, 
                  nalam = integer(1),
                  npass = integer(1), jerr = integer(1), 
                  alsparse = asparse)
  #################################################################################
  # output
  outlist <- getoutput_cpp(fit, maxit, pmax, nvars, vnames)
  if(standardize){
    outlist$beta = outlist$beta/xs
  }
  if(intr){
    outlist$b0 = ym - xm %*% outlist$beta
  }
  outlist <- c(outlist, list(npasses = fit$npass, jerr = fit$jerr, group = group))
  class(outlist) <- c("ls")
  outlist
}

gglasso_cpp <- function(x, y, group = NULL, nlambda = 100, lambda.factor = ifelse(nobs < nvars, 0.01, 1e-04), 
                    lambda = NULL, pf = sqrt(bs), weight = NULL, dfmax = as.integer(max(group)) + 
                      1, pmax = min(dfmax * 1.2, as.integer(max(group))), eps = 1e-08, maxit = 3e+08, 
                    delta, intercept=TRUE, asparse = 0.05, standardize=TRUE) {
  #################################################################################
  #\tDesign matrix setup, error checking
  this.call <- match.call()
  
  if (!is.matrix(x)) 
    stop("x has to be a matrix")
  
  if (any(is.na(x))) 
    stop("Missing values in x not allowed!")
  
  y <- drop(y)
  np <- dim(x)
  nobs <- as.integer(np[1])
  nvars <- as.integer(np[2])
  vnames <- colnames(x)
  
  if (is.null(vnames)) 
    vnames <- paste("V", seq(nvars), sep = "")
  
  if (length(y) != nobs) 
    stop("x and y have different number of rows")
  
  if (!is.numeric(y)) 
    stop("The response y must be numeric. Factors must be converted to numeric")
  
  #################################################################################
  #    group setup
  if (is.null(group)) {
    group <- 1:nvars
  } else if (length(group) != nvars) 
    stop("group length does not match the number of predictors in x")
  
  bn <- as.integer(max(group))
  bs <- as.integer(as.numeric(table(group)))
  
  if (!identical(as.integer(sort(unique(group))), as.integer(1:bn))) 
    stop("Groups must be consecutively numbered 1,2,3,...")
  
  if (asparse>1 || asparse<0){
    asparse = .5
    warning("asparse must be in [0,1], arbitrarily using 0.5.")
  } 
  
  ix <- rep(NA, bn)
  iy <- rep(NA, bn)
  j <- 1
  for (g in 1:bn) {
    ix[g] <- j
    iy[g] <- j + bs[g] - 1
    j <- j + bs[g]
  }
  ix <- as.integer(ix) - 1 # zero indexing
  iy <- as.integer(iy) - 1
  group <- as.integer(group)
  #################################################################################
  #parameter setup
  if (missing(delta)) 
    delta <- 1
  if (delta < 0) 
    stop("delta must be non-negtive")
  delta <- as.double(delta)
  if (length(pf) != bn) 
    stop("The size of group-lasso penalty factor must be same as the number of groups")
  maxit <- as.integer(maxit)
  pf <- as.double(pf)
  eps <- as.double(eps)
  dfmax <- as.integer(dfmax)
  pmax <- as.integer(pmax)
  #################################################################################
  #lambda setup
  nlam <- as.integer(nlambda)
  if (is.null(lambda)) {
    if (lambda.factor >= 1) 
      stop("lambda.factor should be less than 1")
    flmin <- as.double(lambda.factor)
    ulam <- double(1)
  } else {
    #flmin=1 if user define lambda
    flmin <- as.double(1)
    if (any(lambda < 0)) 
      stop("lambdas should be non-negative")
    ulam <- as.double(rev(sort(lambda)))
    nlam <- as.integer(length(lambda))
  }
  intr <- as.integer(intercept)
  #################################################################################
  # call R sub-functions
  fit <- ls_sparse_cpp(bn, bs, ix, iy, nobs, nvars, x, y, pf, 
                   dfmax, pmax, nlam, flmin, ulam, eps, maxit, vnames, group, intr, 
                   asparse, standardize) 
  #################################################################################
  # output
  if (is.null(lambda)) fit$lambda <- gglasso:::lamfix(fit$lambda)
  fit$call <- this.call
  class(fit) <- c("gglasso", class(fit))
  fit
} 

getoutput_cpp <- function(fit, maxit, pmax, nvars, vnames) {
  nalam <- fit$nalam
  nbeta <- fit$nbeta[seq(nalam)]
  nbetamax <- max(nbeta)
  lam <- fit$alam[seq(nalam)]
  stepnames <- paste("s", seq(nalam) - 1, sep = "")
  errmsg <- gglasso:::err(fit$jerr, maxit, pmax)  ### error messages from fortran
  switch(paste(errmsg$n), `1` = stop(errmsg$msg, call. = FALSE), 
         `-1` = print(errmsg$msg, call. = FALSE))
  dd <- c(nvars, nalam)
  if (nbetamax > 0) {
    beta <- fit$beta
    #dimnames(beta) = list(vnames, stepnames)
    #beta <- matrix(fit$beta[seq(nvars * nalam)], nvars, nalam, 
    #dimnames = list(vnames, stepnames))
    df <- apply(abs(beta) > 0, 2, sum)
  } else {
    beta <- Matrix::Matrix(0, nvars, nalam)
    df <- rep(0, nalam)
  }
  b0 <- fit$b0
  if (!is.null(b0)) {
    b0 <- b0[seq(nalam)]
    names(b0) <- stepnames
  }
  list(b0 = b0, beta = beta, df = df, dim = dd, lambda = lam)
}



library(gglasso)
data(bardet)
group1 <- rep(1:20,each=5)
m1 <- gglasso(x=bardet$x,y=bardet$y,group=group1,loss="ls_sparse",asparse = .1)
m2 <- gglasso_cpp(x=bardet$x,y=bardet$y,group=group1,asparse = .1)

*/
