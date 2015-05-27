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


