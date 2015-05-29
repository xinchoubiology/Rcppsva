#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <math.h>
#include <omp.h>

#define OMP_NUM_THREADS omp_get_max_threads()

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
  arma::mat Q, R;
  arma::qr(Q, R, svs);
  Q = Q.cols(0, svs.n_cols-1);
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
    int df = M.n_cols - 1 - svs.n_cols;  // since in QR the rank are all full rank == col(svs)
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


