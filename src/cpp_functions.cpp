// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
double dot_integrate(NumericVector v) {
  double val = 0;

  for(int i=0; i<v.size()-1; i++){
    val += v[i+1]+v[i];
  }
  return val / (2 * (v.size() - 1));
}

// [[Rcpp::export]]
NumericVector dot_integrate_col(NumericMatrix v) {
  NumericVector val(v.ncol());

  for(int i=0; i<v.ncol(); i++){
    val(i) = 0; //v[i][0];
    for(int j=0; j<v.nrow()-1; j++){
      //val[i] += v[i][j+1]+v[i][j];
      val(i) += v(j+1,i)+v(j,i);
    }
    val(i) = val(i) / (2 * (v.nrow() - 1));
  }
  return val;
}

// [[Rcpp::export]]
ComplexMatrix dot_col_cumsum(ComplexMatrix m) {
  for (int i = 1; i < m.nrow(); ++i) {
    for (int j = 0; j < m.ncol(); ++j) {
      m(i, j) = m(i, j) + m(i - 1, j);
    }
  }
  return m;
}

// // [[Rcpp::export]]
// Eigen::MatrixXd eigenMapMatMult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
//   return A * B;
// }
