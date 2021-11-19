// Copyright 2018 Erik-Jan van Kesteren
//   
//   Permission is hereby granted, free of charge, to any person obtaining a 
//   copy of this software and associated documentation files (the "Software"), 
//   to deal in the Software without restriction, including without limitation 
//   the rights to use, copy, modify, merge, publish, distribute, sublicense, 
//   and/or sell copies of the Software, and to permit persons to whom the 
//   Software is furnished to do so, subject to the following conditions:
//     
//   The above copyright notice and this permission notice shall be included 
//   in all copies or substantial portions of the Software.
//   
//   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
//   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
//   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 
//   THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
//   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
//   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
//   DEALINGS IN THE SOFTWARE.


// Resources used:
// http://gallery.rcpp.org/articles/dmvnorm_arma/
// http://gallery.rcpp.org/articles/using-rcppprogress/

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]] 
// [[Rcpp::plugins(openmp)]] 
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <progress.hpp>
#include <progress_bar.hpp>

// Generate random start
arma::uvec generateStart(int n, int p) {
  int nones = std::floor((std::sqrt((double)n) <= p / 2.0) ? std::sqrt((double)n) : p / 2.0);
  arma::uvec one = arma::uvec(nones, arma::fill::ones);
  arma::uvec zer = arma::uvec(p - nones, arma::fill::zeros);
  return arma::shuffle(arma::join_cols(one, zer));
}


// Column filter function
arma::mat getcols(arma::mat & M, arma::uvec & idx) {
  int p = arma::sum(idx);
  arma::mat A = arma::mat(M.n_rows, p);
  int index = 0;
  for (arma::uword i = 0; i < M.n_cols; i++) {
    if (idx(i) == 1) {
      arma::vec candidate = M.col(i);
      A.col(index) = candidate;
      index++;
    }
  }
  return A;
}


// fast implementation of sobel test
double sobel(arma::vec & x, arma::vec & m, arma::vec & y) {
  // get alpha and its variance
  int n           = x.n_elem;
  double cpx      = dot(x,x);
  double alpha    = dot(x,m) / cpx;
  arma::vec res_m = m - x*alpha;
  double var_a    = dot(res_m, res_m) / (n - 1) / cpx;
  
  // construct model matrix for b
  arma::mat mm    = arma::mat(n, 2);
  mm.col(0)       = x;
  mm.col(1)       = m;
  arma::mat cpmm  = mm.t() * mm;
  arma::vec cpmy  = mm.t() * y;
  
  // get beta and its variance
  arma::mat cpmmi = cpmm.i();
  arma::vec beta  = cpmmi * cpmy;
  arma::vec res_y = y - mm * beta;
  arma::mat var_b = dot(res_y, res_y) / (n - 2) * cpmmi;
  
  // return the z-score
  double stat = alpha * beta(1);
  double se   = std::sqrt(alpha*alpha*var_b(1,1) + beta(1)*beta(1)*var_a);
  return stat/se;
}

// fast implementation of causal steps test
double csteps(arma::vec & x, arma::vec & m, arma::vec & y) {
  // get alpha and its variance and t score
  int n           = x.n_elem;
  double cpx      = dot(x,x);
  double alpha    = dot(x,m) / cpx;
  arma::vec res_m = m - x*alpha;
  double var_a    = dot(res_m, res_m) / (n - 1) / cpx;
  double ta       = fabs(alpha / std::sqrt(var_a));
  
  // construct model matrix for b
  arma::mat mm    = arma::mat(n, 2);
  mm.col(0)       = x;
  mm.col(1)       = m;
  arma::mat cpmm  = mm.t() * mm;
  
  // get beta and its variance and t score
  arma::vec beta  = arma::solve(cpmm, mm.t() * y);
  arma::vec res_y = y - mm * beta;
  arma::mat var_b = dot(res_y, res_y) / (n-1) * cpmm.i();
  double tb       = fabs(beta(1) / std::sqrt(var_b(1, 1)));
  
  // return min(a, b)
  return (ta < tb) ? ta : tb;
}

// fast implementation of cmf step
arma::uvec cmfastep(arma::vec & x, arma::mat & M, arma::vec & y, 
                    arma::uvec f, arma::uvec & s, 
                    double & critval, int & decfun) {
  // Rcpp::Rcout << 1 << std::endl;
  int n = x.n_elem;                                // sample size
  int sl = s.n_elem;                               // n Ms to consider
  arma::vec m;                                     // init the mediator
  arma::vec r_x;                                   // init residual of x
  arma::vec r_y;                                   // init residual of y

  for (int i = 0; i < sl; i++) {
    int idx = s(i);                                // idx of M to consider
    int curcsum = arma::sum(f) - f(idx);           // n remaining selected M
    if (curcsum >= n-2) {
      Rcpp::Rcout << "c";
      break;
    }                  // escape if p > n
    m = M.col(idx);                                // assign M[, idx]
    if (curcsum == 0) {                            // no remaining Ms selected
      r_x = x;
      r_y = y;
    } else {
      arma::uvec fcopy = f;                        // create copy of filter
      fcopy(idx)       = 0;                        // without current M
      
      arma::mat Mx   = getcols(M, fcopy);          // model matrix
      arma::mat MtM  = Mx.t() * Mx;                // crossprod
      //arma::mat icp  = arma::inv_sympd(MtM);        // inverse crossprod matrix
      arma::mat H    = Mx * arma::solve(MtM, Mx.t(), arma::solve_opts::likely_sympd); // hat matrix
      r_x            = x - H * x;                  // residual of x
      r_y            = y - H * y;                  // residual of y
    }
    
    // test for 
    bool isMediator;
    if (decfun == 1) {
      isMediator = fabs(sobel(r_x, m, r_y)) > critval; // perform sobel test
    } else {
      isMediator = csteps(r_x, m, r_y) > critval;      // causal steps test
    }
    
    if (isMediator) {
      f(idx) = 1;                                 // select m
    } else {
      f(idx) = 0;                                 // deselect m
    }
  }
  
  return f;
}


// fast implementation of cmf
// [[Rcpp::export]]
arma::vec arma_cmf(arma::vec & x, arma::mat & M, arma::vec & y, 
                   int & maxIter, int & stableLag, 
                   double & critval, int & decfun,
                   int & nCores, int & nStarts, bool & pb) {
  // set the number of cores to use
#ifdef _OPENMP
  omp_set_num_threads(nCores);
#endif
  
  // initialise information
  int n = M.n_rows;
  int p = M.n_cols;
  int ss = std::ceil(std::sqrt((double)p));
  
  // initialise output matrix to fill with selections
  arma::umat result = arma::umat(M.n_cols, nStarts, arma::fill::zeros);
  
  // instantiate progress bar
  Progress prog(nStarts, pb);
  
  // start multithreading
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (int start = 0; start < nStarts; start++) {
    // generate random start for the filter vector
    arma::uvec f = generateStart(n, p);
    if ( !Progress::check_abort() ) {
      if (pb) prog.increment(); // update progress
      
      // prepare for inner loop
      int i = 0; // index
      arma::uvec subset(p); // subsetting vector
      
      // prepare convergence matrix 
      bool converged = false;
      arma::umat conv(p, stableLag + 1, arma::fill::ones);
      conv.col(stableLag) = f;
      
      // inner loop
      while ((i < maxIter) & !converged) {
        // get subset in random order
        subset = arma::shuffle(arma::linspace<arma::uvec>(0, p-1, p)); 
        arma::uvec s = subset.head(ss);
        
        // update the filter vector
        f = cmfastep(x, M, y, f, s, critval, decfun); 
        
        // update convergence matrix
        conv.shed_col(0);
        conv = arma::join_rows(conv, f);
        
        // check for convergence
        int ccount = p;
        for (int j = 0; j < p; j++) {
          int sumcrow = sum(conv.row(j));
          ccount -= sumcrow % (stableLag + 1);
        }
        if (ccount == p) converged = true;
        i++;
      }
    }
    
    result.col(start) = f;
  }
  
  // reduce to empirical selection probabilities and return
  return arma::conv_to<arma::vec>::from(arma::sum(result, 1)) / nStarts;
}

// Check if openMP is available
// [[Rcpp::export]]
bool checkOMP() {
#ifdef _OPENMP
  return true;
#else
  return false;
#endif
}
