#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <math.h>
#include <vector>
#include <string>
#include <omp.h>

#define OMP_NUM_THREADS omp_get_max_threads()

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]


//' Integrating empirical bayesian prior function
//' by Monte Carlo integration
//'
//' @title int_eprior
//' @param sx standard data matrix
//' @param ghat estimated gamma values
//' @param dhat estimated delta values
//' @return list
//'         g.star
//'         d.star
//' @author Xin Zhou \url{xinchoubiology@@gmail.com}
// [[Rcpp::export]]
Rcpp::List int_eprior(const Eigen::MatrixXd & sx, const Eigen::VectorXd & ghat, const Eigen::VectorXd & dhat){
  int r = sx.rows(), c = sx.cols();

  Eigen::VectorXd gstar = Eigen::VectorXd::Zero(r);
  Eigen::VectorXd dstar = Eigen::VectorXd::Zero(r);

  #pragma omp parallel for schedule(dynamic, 32) num_threads(OMP_NUM_THREADS)
  for(int i = 0; i < r; i++){
    Eigen::MatrixXd dat  = Eigen::MatrixXd::Zero(r, c);
    dat.rowwise()        = sx.row(i);
    dat.colwise()       -= ghat;
    dat                  = dat.array().square();
    Eigen::MatrixXd sum2 = dat.rowwise().sum();     // dat %*% I
    Eigen::MatrixXd LH   = 1.0 /  (sqrt(pow(2 * M_PI, c)) * dhat.array().pow(c/2.0)) * \
                           (-1 * (sum2.array() / (2 * dhat.array()))).exp();
    Eigen::MatrixXd gLH  = LH.array() * ghat.array();
    Eigen::MatrixXd dLH  = LH.array() * dhat.array();
    gstar(i)             = (gLH.sum() - gLH(i)) / (LH.sum() - LH(i));
    dstar(i)             = (dLH.sum() - dLH(i)) / (LH.sum() - LH(i));
  }
  return Rcpp::List::create(Rcpp::Named("g.star") = gstar,
                            Rcpp::Named("d.star") = dstar);
}

//' SVD calculation for X=UΣV, use XtX=VtΣV to calculate the right eigen vector
//' @title arma_eigen
//' @param M input matrix = XtX
//' @return v eigen vector
//' @export
//' @author Xin Zhou \url{xinchoubiology@@gmail.com}
//  [[Rcpp::export]]
arma::vec arma_eigen(const arma::mat & M){
  return arma::eig_sym(M);
}


//' Linear regression on covariate of interest when comparing with other covariates
//' @title beta_regress
//' @param M m x n expression matrix; Each row represents probes and each col means samples
//' @param pv n x B design matrix; Each col means phenotype of interest, and if B >=2, 
//'        means we have B-1 permutations on covariates of interest
//' @param svs n x (p-1) design matrix; Each row represent samples and each col means parameters
//' @param full full output or coefficient only
//' @export
//' @author Xin Zhou
// [[Rcpp::export]]
Rcpp::List beta_regress(const arma::mat & M, const arma::mat & pv, const arma::mat & svs, const int full = 0){
  int n = svs.n_rows;
  int m = svs.n_cols;
  arma::mat Q, R;
  arma::qr(Q, R, svs);
  Q = Q.cols(0, arma::rank(R) - 1);                           // svs.n_cols-1 is Q's rank
  arma::mat S   = arma::eye<arma::mat>(n, n) - Q * Q.t();     // n x n
  arma::mat sv  = S * pv;                                     // sv n x n x n x B = n x B
  arma::vec vsv = arma::diagvec(arma::trans(pv) * sv);        // B x n x n x B = B x B
  arma::mat b   = M * arma::trans(S) * pv;                    // m x n x n x n x n x B = m x B
  b.each_row() /= vsv.t();
  if(full == 0){
    return Rcpp::List::create(Rcpp::Named("coef") = b);
  }
  else{
    arma::mat sy = M * S;    // m x n x n x n = m x n
    int df = M.n_cols - 1 - m;  // since in QR the rank are all full rank == col(svs)
    int B = pv.n_cols;
    arma::mat sigma(M.n_rows, B);
    if(B == 1){
      sigma = arma::sqrt(arma::sum(arma::square(sy - b * sv.t()), 1) / df);
    }
    else{
      #pragma omp parallel for schedule(dynamic, 32) shared(sy) num_threads(OMP_NUM_THREADS)
      for(int i = 0; i < B; i++){
        sigma.col(i) = arma::sum(arma::square(sy - b.col(i) * sv.col(i).t()), 1);
      }
      sigma = arma::sqrt(sigma / df);
    }
    return Rcpp::List::create(Rcpp::Named("coef") = b,
                              Rcpp::Named("sigma") = sigma,
                              Rcpp::Named("stdev_unscaled") = arma::sqrt(1 / vsv),
                              Rcpp::Named("df.residuals") = df);
  }
}

//' Bootstrap testing on regression model for covariate of interested
//' 
//' @title bootstrap_regress
//' @description Dimitris[Bootstrap hypothesis testing in regression models] 
//' @param M m x n expression matrix; Each row represents probes and each col means samples
//' @param mod n x p design matrix; 
//' @param modn n x (p-1) null design matrix - covariate of interest;
//' @param B n x B matrix; Bootstrap iterations index matrix
//' @export
//' @author Xin Zhou
//  [[Rcpp::export]]
Rcpp::List bootstrap_regress(const arma::mat & M, const arma::mat & mod, const arma::mat & modn, const arma::umat & B){
  int n = mod.n_rows;
  arma::mat Q, R;
  arma::qr(Q, R, mod);
  int rank = arma::rank(R);
  Q = Q.cols(0, rank - 1);
  R = R.rows(0, rank - 1);
  arma::mat S   = arma::eye<arma::mat>(n, n) - Q * Q.t();     // n x n
  arma::mat resid = M * S;                                    // m x n x n x n = m x n
  arma::vec scale = arma::sqrt(arma::diagvec(S));             // n x 1
  resid.each_row() /= scale.t();
  // calculate null M
  arma::mat Qn, Rn;
  arma::qr(Qn, Rn, modn);
  Qn = Qn.cols(0, arma::rank(Rn) - 1);
  arma::mat null = M * Qn * Qn.t();                          // m x n x n x n = m x n
  // Bootstrap test
  arma::mat beta0  = arma::zeros<arma::mat>(M.n_rows, B.n_cols);    // m x B
  arma::mat sigma0 = arma::zeros<arma::mat>(M.n_rows, B.n_cols);    // m x B
  int df = M.n_cols - mod.n_cols;
  
  #pragma omp parallel for schedule(dynamic, 32) num_threads(OMP_NUM_THREADS)
  for(unsigned int i = 0; i < B.n_cols; i++){
    arma::uvec index = arma::vectorise(B.col(i));             // n vector
    arma::mat Mb  = null + resid.cols(index);                 // m x n + m x n = m x n
    arma::mat tmp = arma::inv(R) * Q.t() * Mb.t();            // p x p x p x n x n x m = p x m
    beta0.col(i)  = tmp.row(1).t();                           // m vector assignment
    // m x 1
    sigma0.col(i) = arma::sum(arma::square(Mb - Mb * Q * Q.t()), 1); 
  }
  sigma0 = arma::sqrt(sigma0 / df);
  
  return Rcpp::List::create(Rcpp::Named("coef") = beta0,
                            Rcpp::Named("sigma") = sigma0,
                            Rcpp::Named("df.residuals") = df);
}


//' Parallel Computing measurement matrix for different linkage type : single linkage
//' @title single_linkage
//' @param M correlation matrix
//' @return Matrix
//' @export
//  [[Rcpp::export]]
arma::mat single_linkage(arma::mat & M){
  arma::mat N = M;
  for(unsigned i = 0; i < N.n_rows; i++){
    for(unsigned j = i+2; j < N.n_cols; j++){
      N(i, j) = Rcpp::min(Rcpp::NumericVector::create(N(i, j-1), N(i+1,j), N(i,j)));
    }
  }
  return N;
}


//' Parallel Computing measurement matrix for different linkage type : complete linkage
//' @title complete_linkage
//' @param M correlation matrix
//' @return Matrix
//' @export
//  [[Rcpp::export]]
arma::mat complete_linkage(arma::mat & M){
  arma::mat N = M;
  for(unsigned i = 0; i < N.n_rows; i++){
    for(unsigned j = i+2; j < N.n_cols; j++){
      N(i, j) = Rcpp::max(Rcpp::NumericVector::create(N(i, j-1), N(i+1,j), N(i,j)));
    }
  }
  return N;
}


//' Parallel Computing measurement matrix for different linkage type : average linkages
//' @title average_linkage
//' @param M correlation matrix
//' @return Matrix
//' @export
//  [[Rcpp::export]]
arma::mat average_linkage(arma::mat & M){
  arma::mat N = M;
  for(unsigned i = 0; i < N.n_rows; i++){
    for(unsigned j = i+2; j < N.n_cols; j++){
      N(i, j) += N(i, j-1);
      for(unsigned k = i+1; k < j; k++){
        N(i, j) += N(k, j);
      }
    }
  }
  
  #pragma omp parallel for num_threads(OMP_NUM_THREADS)
  for(unsigned i = 0; i < N.n_rows; i++){
    for(unsigned j = i+2; j < M.n_cols; j++){
      N(i, j) /= ((j - i) * (j - i + 1) / 2);
    }
  }
  return N;
}

//' Parallel clique merge by indicator matrix by prefix-product
//' 
//' @title clique_merge
//' @param M Indicator Square Matrix
//' @return index NumericVector
//' @export
//  [[Rcpp::export]]
Rcpp::NumericVector clique_merge(arma::mat M){
  int N = M.n_cols;
  arma::mat T(M.n_rows, M.n_cols);
  arma::mat K(M.n_rows, M.n_cols);
  Rcpp::NumericVector res(N - 1);
  for(int n = N - 1; n > 1; n--){
    for(int j = 0; j < std::log2(n - 1); j++){
      #pragma omp parallel for num_threads(OMP_NUM_THREADS)
      for(int i = 1<<j; i < n; i++){
        T(N - 1 - n, n - i) = M(N - 1 - n, n - i) * M(N - 1 - n, n - i + 1);
        K(i, n) = M(i, n) * M(i - 1, n);
      }
      
      #pragma omp parallel for num_threads(OMP_NUM_THREADS)
      for(int i = 1<<j; i < n; i++){
        M(N - 1 - n, n - i) = T(N - 1 - n, n - i);
        M(i, n) = K(i, n);
      }
    }
    #pragma omp parallel for num_threads(OMP_NUM_THREADS)
    for(int k = N - n + 1; k < n; k++){
      M(N - n, k) *= M(N - 1 - n, k);
      M(N - k, n - 1) *= M(N - k, n);
    }
  }

  #pragma omp parallel for num_threads(OMP_NUM_THREADS)
  for(int i = 0; i < N - 1; i++){
    res[i] = M(i, i + 1);
  }
  
  return res;
}
