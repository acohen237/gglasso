
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;


void softthresh(arma::vec x, double thresh){
  arma::vec::iterator it = x.begin();
  arma::vec::iterator it_end = x.end();
  int sgn;
  double shrink;
  for(;it != it_end; ++it){
    sgn = arma::sign(*it);
    shrink = abs(*it) - thresh;
    shrink = (shrink > 0) ? shrink : 0;
    *it = sgn*shrink;
  }
}

void ls_f_sparse(int bn, arma::uvec bs, arma::uvec ix,arma::uvec iy,arma::vec gam,
                 int nobs, int nvars, arma::mat x, arma::vec y,
                 arma::vec pf, int dfmax, int pmax, int nlam, double flmin,
                 arma::vec ulam, double eps, int maxit, int intr, int nalam,
                 arma::vec b0, arma::mat beta, arma::uvec idx,
                 arma::uvec nbeta, arma::vec alam, int npass, int jerr,
                 double alsparse){
  double big=9.9E30;
  double mfl = 1.0E-6;
  int mnlam = 6;
  int mnl;
  
  // - - - local declarations - - -
  double max_gam, d, dif, al, alf, sg;
  arma::vec b(nvars);
  arma::vec oldbeta(nvars);
  arma::vec r = y; // Residue y-beta_k*x etc
  arma::vec oldb;
  arma::vec dd;
  arma::vec oidx(bn);
  int g, j, l, ni, me, startix, endix;
  // - - - Aaron's declarations
  double snorm;
  arma::vec t_for_s(bn); // this is for now just 1/gamma
  double tea; // this takes the place of 't' in the update step for ls
  arma::vec s; // takes the place of 'u' in update for ls
  int vl_iter; // for iterating over columns(?) of x*r
  int kill_count=0;
               
  double tlam, lama, lam1ma, al0=0.0, testr=0.0, tmp;
  int jx;
  arma::uvec jxx(bn);
  arma::vec ga(bn);
  arma::vec vl(nvars);
  max_gam = gam.max();
  if(max(pf) <= 0.0){
    jerr=10000;
    return;
  }
  pf = arma::clamp(pf, 0.0, pf.max());
  jxx.zeros();
  al = 0.0;
  mnl = (nlam >= mnlam) ? mnlam : nlam; // min
  b.zeros();
  oldbeta.zeros();
  idx.zeros();
  oidx.zeros();
  npass = 0;
  ni = npass;
  alf = 0.0;
  t_for_s = 1/gam;
  // --------- lambda loop ----------------------------
  if(flmin < 1.0){ // THIS is the default...
    flmin = (mfl >= flmin) ? mfl : flmin; // just sets a threshold above zero 
    alf = pow(flmin, 1.0/(nlam-1.0));
  }
  // PRINT *, alf
  vl = x.t() * r / nobs; // Note r gets updated in middle and inner loop 
  al0 = 0.0;
  // For each group...      
  for(g = 0; g<bn; g++) ga(g) = arma::norm(vl.subvec(ix(g),iy(g)));
  al0 = arma::norm(vl, "inf"); // Infty norm of X'y, big overkill for lam_max
  // PRINT *, alsparse
  al = al0; // ! this value ensures all betas are 0
  l = 0;
  tlam = 0.0;
  while (l < nlam){ // This is the start of the loop over all lambda values...
    // print *, "l = ", l
    // print *, "al = ", al
    // IF(kill_count > 1000) RETURN
    // kill_count = kill_count + 1
    al0 = al; // store old al value on subsequent loops, first set to al
    if (flmin>=1.0){ // user supplied lambda value, break out of everything
      l++;
      al=ulam(l);
      // print *, "This is at the flmin step of the while loop"
    } else {
      if (l > 1){ // have some active groups
	      al*=alf;
        tmp = 2.0*al-al0;
	      tlam = (tmp >= 0) ? tmp : 0.0; // Here is the strong rule...
	      l++;
	      // print *, "This is the l>1 step of while loop"
      } else if (l==0){ //Trying to find an active group
	      al*=.99;
	      tlam=al;
      }
    }
    lama = al*alsparse;
    lam1ma = al*(1-alsparse);
    // This is the start of the algorithm, for a given lambda...
    // print *, "Here is tlam = ", tlam
    for(g = 0; g < bn; g++){ 
      if(jxx(g) == 1)  continue;
      if(ga(g) > pf(g)*tlam*(1-alsparse)) jxx(g)++; // Implementing the strong rule
    }
    // ! --------- outer loop ---------------------------- ! 
    do { // while condition checked way down below (jx==1)
      oldbeta(0)=b(0);
      // print *, "Here is the outer loop, and here's oldbeta:", oldbeta
      if(ni>0){
	      for(j=0; j<ni; j++){
	        g=idx(j);
	        oldbeta.subvec(ix(g),iy(g))=b.subvec(ix(g),iy(g));
	      }
      }
      // --middle loop-------------------------------------
      for(;;) { // loops forever, we rely on one of the break statements
	      // print *, "This is where we enter the middle loop"
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
	        if(!dd.is_zero()){
	          tmp = gam(g)*gam(g)*arma::dot(dd,dd);
	          dif = (dif >= tmp) ? dif : tmp;
	          r -= x.cols(startix, endix) * dd;
	          if(oidx(g)==0){// Here is where middle loop is different;
	            // if group g was not in oidx (active), and the
	            // difference was nonzero, put it in active (ni)
	            ni++;
	            if(ni>pmax) break;
	            oidx(g)=ni;
	            idx(ni)=g;
	          }
	        }
	      }
	      if(intr != 0){
	        d = arma::mean(r);
	        if(d != 0.0){
	          b(0)+=d;
	          r -= d;
	          tmp = d*d;
	          dif = (dif >= tmp) ? dif : tmp;
	        }
	      }
	      if (ni > pmax) break;
	      if (dif < eps) break;
	      if (npass > maxit){
	        jerr = -l;
	        return;
	      }
	      // --inner loop----------------------
	      do { // do...while() loop, see the end for the condition
	        // PRINT *, "Here is where the inner loop starts"
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
	        if(intr != 0){
	          d = arma::mean(r);
	          if(d != 0.0){
	            b(0)+=d;
	            r -= d;
	            tmp = d*d;
	            dif = (dif > tmp) ? dif : tmp;
	          }
	        }
	        if(npass > maxit){
	          jerr = -l;
	          return;
	        }
	      } while(dif < eps);  // End inner loop
      } // End middle loop                    
      if(ni>pmax) break;
      // !--- final check ---------- ! This checks which violate KKT condition
      // PRINT *, "Here is where the final check starts"
      jx = 0;
      for(j=0; j<nvars; j++){
        testr = max_gam*(b(j)-oldbeta(j))/(1+abs(b(j)));
        testr *= testr;
        if(testr >= eps){
          jx = 1;
          break;
        }
      }
      if (jx == 1) continue; // return to top of outer loop
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
    } while(jx==1); // Ends outer loop
    // ---------- final update variable and save results------------
    if (l==0){
      if(jxx.is_zero()){
        continue; // don't save anything, we're still decrementing lambda
      } else {
        l=2;
        tmp = (alf >= .99) ? alf : .99;
        alam(0) = al / tmp; // store previous, larger value
      }
    }
    // PRINT *, "Here is where the final update starts"
    if (ni>pmax) {
      jerr = -10000-l;
      break;
    }
    if (ni>0) {
      for(j=0; j<ni; j++){
        g=idx(j);
        beta(arma::span(ix(g),iy(g)),l-1)=b.subvec(ix(g),iy(g));
      }
    }
    nbeta(l-1)=ni;
    b0(l-1)=b(0);
    alam(l-1)=al;
    nalam=l;
    if (l < mnl) continue;
    me=0;
    for(j=0; j<ni; j++){
      g=idx(j);
      if(! beta(arma::span(ix(g),iy(g)),l-1).is_zero()) me++;
    }
    if (me>dfmax) break;
  } // end lambda loop
  return;
}
