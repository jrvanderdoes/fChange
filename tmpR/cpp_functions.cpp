#include "cpp11.hpp"
using namespace cpp11;

[[cpp11::register]]
double dot_integrate(doubles v) {
    double val = 0;

    for(int i=0; i<v.size()-1; i++){
        val += v[i+1]+v[i];
    }
    return val / (2 * (v.size() - 1));
}

[[cpp11::register]]
doubles dot_integrate_col(doubles_matrix<> v) {
  writable::doubles val(v.ncol());

  for(int i=0; i<v.ncol(); i++){
    val[i] = 0; //v[i][0];
    for(int j=0; j<v.nrow()-1; j++){
      val[i] += v[i][j+1]+v[i][j];
    }
    val[i] = val[i] / (2 * (v.nrow() - 1));
  }
  return val;
}
