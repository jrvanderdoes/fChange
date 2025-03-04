// [[Rcpp::depends(RcppArmadillo)]]
#include "RcppArmadillo.h"

using namespace Rcpp;

// // [[Rcpp::export]]
// double dot_integrate(NumericVector v) {
//   double val = 0;
//
//   for(int i=0; i<v.size()-1; i++){
//     val += v[i+1]+v[i];
//   }
//   return val / (2 * (v.size() - 1));
// }

// [[Rcpp::export]]
double dot_integrate(NumericVector v,
                     Nullable<NumericVector> r = R_NilValue) {
  double val = 0;
  NumericVector rr;

  // Make sure we have the proper resolution
  if(r.isUsable()) {

    rr = r;

  } else{

    rr = NumericVector(v.size());
    for (int i = 0; i < rr.length(); i++) {
      rr(i) = i * ( 1.0 / (rr.length()-1) );
    }

  }

  for(int i=0; i<v.size()-1; i++){
    val += ( v(i+1)+v(i) ) * ( rr(i+1)-rr(i) );
  }
  return val / 2 ;
}

// // [[Rcpp::export]]
// NumericVector dot_integrate_col(NumericMatrix v) {
//   NumericVector val(v.ncol());
//
//   for(int i=0; i<v.ncol(); i++){
//     val(i) = 0; //v[i][0];
//     for(int j=0; j<v.nrow()-1; j++){
//       //val[i] += v[i][j+1]+v[i][j];
//       val(i) += v(j+1,i)+v(j,i);
//     }
//     val(i) = val(i) / (2 * (v.nrow() - 1));
//   }
//   return val;
// }

// [[Rcpp::export]]
NumericVector dot_integrate_col(NumericMatrix v,
                                Nullable<NumericVector> r = R_NilValue) {
  NumericVector val(v.ncol());
  NumericVector rr;

  // Make sure we have the proper resolution
  if(r.isUsable()) {

    rr = r;

  } else{

    rr = NumericVector(v.nrow());
    for (int i = 0; i < rr.length(); i++) {
      rr(i) = i * ( 1.0 / (rr.length()-1) );
    }

  }

  // Integrate
  for(int i=0; i<v.ncol(); i++){
    val(i) = 0;
    for(int j=0; j<v.nrow()-1; j++){ // resolution
      val(i) += ( v(j+1,i)+v(j,i) ) * ( rr(j+1)-rr(j) );
    }
    val(i) = val(i) / 2;
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

// [[Rcpp::export]]
arma::mat dot_sqrt_mat(arma::mat A){
//arma::cx_mat dot_sqrt_mat(arma::mat A){

  return arma::sqrtmat_sympd(A);
  //return arma::sqrtmat(A);
}

