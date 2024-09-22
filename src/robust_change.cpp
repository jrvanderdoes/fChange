// [[Rcpp::depends(RcppArmadillo)]]
#include "RcppArmadillo.h"
#include <vector>
#include <thread>

using namespace Rcpp;

NumericMatrix make_Obs_tilde_h(NumericMatrix Obs);
NumericMatrix make_Obs_tilde_C(NumericMatrix Obs);
double adaptive_bw(NumericMatrix Obs);
NumericMatrix matadd(NumericMatrix A, NumericMatrix B);
NumericVector kernel(double bw, int n);
NumericMatrix toeplitz(NumericVector v);
NumericMatrix getRealSQM(NumericMatrix M);
NumericMatrix make_hC_Obs(NumericMatrix Obs);
NumericMatrix fill_U(NumericMatrix A_h, NumericMatrix A_C, NumericVector norm, NumericMatrix hC_Obs, int m, int n, int d);
NumericVector vecmult(NumericMatrix M, NumericVector v);
NumericMatrix outerProd(NumericVector v, NumericVector w);
NumericVector find_rowmax(NumericMatrix U_hC);
NumericMatrix fill_T(NumericMatrix hC_Obs, int n, int d);
NumericVector h_cpp(NumericVector v1, NumericVector v2);


/// @param[in] nb_elements : size of your for loop
/// @param[in] functor(start, end) :
/// your function processing a sub chunk of the for loop.
/// "start" is the first index to process (included) until the index "end"
/// (excluded)
/// @code
///     for(int i = start; i < end; ++i)
///         computation(i);
/// @endcode
/// @param use_threads : enable / disable threads.
///
///
static
void parallel_for(unsigned nb_elements,
                  std::function<void (int start, int end)> functor,
                  bool use_threads = true) {
  // -------
  unsigned nb_threads_hint = std::thread::hardware_concurrency();
  unsigned nb_threads = nb_threads_hint == 0 ? 8 : (nb_threads_hint);

  unsigned batch_size = nb_elements / nb_threads;
  unsigned batch_remainder = nb_elements % nb_threads;

  std::vector< std::thread > my_threads(nb_threads);

  if( use_threads )
  {
    // Multithread execution
    for(unsigned i = 0; i < nb_threads; ++i)
    {
      int start = i * batch_size;
      my_threads[i] = std::thread(functor, start, start+batch_size);
    }
  }
  else
  {
    // Single thread execution (for easy debugging)
    for(unsigned i = 0; i < nb_threads; ++i){
      int start = i * batch_size;
      functor( start, start+batch_size );
    }
  }

  // Deform the elements left
  int start = nb_threads * batch_size;
  functor( start, start+batch_remainder);

  // Wait for the other thread to finish their task
  if( use_threads )
    std::for_each(my_threads.begin(), my_threads.end(), std::mem_fn(&std::thread::join));
}


// for Wilcoxon-type: calculate \tilde(X)_i = 1/(n-1) \sum_{j=1,\neq j}^n h(X_i,X_j) for data adapted bandwidth
// [[Rcpp::export]]
NumericMatrix make_Obs_tilde_h(NumericMatrix Obs){
  double n = Obs.ncol();
  double d = Obs.nrow();
  NumericMatrix Obs_tilde(d,n);
  NumericVector v1(d);
  NumericVector v2(d);

  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      if(j!=i){
        for(int D=0; D<d; D++){
          v1(D)=Obs(D,i);
          v2(D)=Obs(D,j);
        }
        NumericVector help = h_cpp(v1, v2);
        for(int D=0; D<d; D++){
          Obs_tilde(D,i) = Obs_tilde(D,i) + help(D);
        }
      }
    }
  }
  return Obs_tilde/(n-1);
}


// for CUSUM: calculate \tilde(X)_i = 1/(n-1) \sum_{j=1,\neq j}^n (X_i - X_j) for data adapted bandwidth
// [[Rcpp::export]]
NumericMatrix make_Obs_tilde_C(NumericMatrix Obs){
  double n = Obs.ncol();
  double d = Obs.nrow();
  NumericMatrix Obs_tilde(d,n);

  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      if(j!=i){
        for(int D=0;D<d;D++){
          Obs_tilde(D,i)=Obs_tilde(D,i)+(Obs(D,i)-Obs(D,j));
        }
      }
    }
  }
  return Obs_tilde/(n-1);
}


// evaluate data-adaptive bandwidth
// [[Rcpp::export]]
double adaptive_bw(NumericMatrix Obs){
  double n=Obs.ncol();
  double d=Obs.nrow();
  NumericMatrix autoVs(d,d);
  NumericMatrix CP_0(d,d);
  NumericMatrix CP_1(d,d);
  NumericVector v(d);
  NumericVector w(d);

  double bw_0=ceil(pow(n,(1.0/5.0)));
  //double bw_0=ceil(pow(n,(1.0/3.0)));
  for(double k=0; k<bw_0;k++){
    for(int D=0;D<d;D++){
      v(D)=Obs(D,0);
      w(D)=Obs(D,k);
    }
    autoVs=outerProd(v,w);
    for(int i=1;i<(n-k);i++){
      for(int D=0; D<d;D++){
        v(D)=Obs(D,i);
        w(D)=Obs(D,(k+i));
      }
      autoVs=matadd(autoVs,outerProd(v,w));
    }
    autoVs=autoVs/n;
    if(k==0){
      CP_0=autoVs;
    }
    //CP_0=matadd(CP_0,2.0*(1.0-(k/bw_0))*autoVs);  //Bartlett kernel
    //CP_1=matadd(CP_1,2.0*k*(1.0-(k/bw_0))*autoVs);
    else{
      double quadspec = (25/(12*M_PI*M_PI*k*k/(bw_0*bw_0)))*(sin(6*M_PI*(k/bw_0)/5)/(6*M_PI*(k/bw_0)/5)-cos(6*M_PI*(k/bw_0)/5));
      CP_0=matadd(CP_0,2.0*quadspec*autoVs);
      CP_1=matadd(CP_1,2.0*k*quadspec*autoVs);
    }
  }
  double CPdiag=0;
  for(int D=0;D<d;D++){
    CPdiag=CPdiag+CP_0(D,D);
  }
  double sum1=0;
  double sum0=0;

  for(int i=0;i<d;i++){
    for(int j=0;j<d;j++){
      sum1=sum1+(CP_1(i,j)*CP_1(i,j));
      sum0=sum0+(CP_0(i,j)*CP_0(i,j));
    }
  }
  double bw_opt = pow(n*3.0*sum1/(sum0+(CPdiag*CPdiag)),(1.0/5.0));
  return bw_opt=ceil(bw_opt);

}


// [[Rcpp::export]]
NumericMatrix matadd(NumericMatrix A, NumericMatrix B){
  int m=A.ncol();
  int n=A.nrow();
  NumericMatrix M(n,m);

  for(int i=0;i<n;i++){
    for(int j=0;j<m;j++){
      M(i,j)=A(i,j)+B(i,j);
    }
  }
  return M;
}


// [[Rcpp::export]]
NumericVector kernel(double bw, int n){
  NumericVector ker(n);
  ker(0)=1;
  for(double i=1; i<n; i++){
    ker(i)=(25/(12*M_PI*M_PI*i*i/(bw*bw)))*(sin(6*M_PI*(i/bw)/5)/(6*M_PI*(i/bw)/5)-cos(6*M_PI*(i/bw)/5));
  }
  return(ker);
}


// Toeplitz norm of Vector v
// [[Rcpp::export]]
NumericMatrix toeplitz(NumericVector v){
  int n = v.length();
  NumericMatrix KQS(n,n);
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      KQS(i,j)= v((std::abs(i-j)));
    }
  }
  return KQS;
}


// evaluate real part of matrix square root
// [[Rcpp::export]]
NumericMatrix getRealSQM(NumericMatrix M){
  arma::mat M_arma = as<arma::mat>(M);
  arma::cx_mat SQM = arma::sqrtmat(M_arma);
  arma::mat SQM_real_arma = real(SQM);
  NumericMatrix SQM_real = as<Rcpp::NumericMatrix>(wrap(SQM_real_arma));
  return SQM_real;
}


// calculate h(X_i,X_j) and X_i - X_j and save in one matrix
// [[Rcpp::export]]
NumericMatrix make_hC_Obs(NumericMatrix Obs){
  int d = Obs.nrow();
  int n = Obs.ncol();
  NumericMatrix hC_Obs(d*(n-1),n);//2*n);
  NumericVector v1(d);
  NumericVector v2(d);


  for (int i=0; i<n-1; i++){
    for (int j=i+1; j<n; j++){
      for (int D=0; D<d; D++){
        v1(D)=Obs(D,i);
        v2(D)=Obs(D,j);
      }
      NumericVector help = h_cpp(v1, v2);
      for (int D=0; D<d; D++){
        hC_Obs(D*(n-1)+i,j) = help(D);  //entries for h_Obs
        // hC_Obs(D*(n-1)+i,j+n)= Obs(D,i)-Obs(D,j); //entries for C_Obs
      }
    }
  }

  return hC_Obs;
}


// all of the above calculations are done once for simulated observations X_1,...,X_n.
// for each bootstrap iteration t=0,...,m, add multiplier and calculate U_{n,k}^(t) for k=1,...,n
// [[Rcpp::export]]
NumericMatrix fill_U(NumericMatrix A_h, //NumericMatrix A_C,
                     NumericVector norm, NumericMatrix hC_Obs,
                     int m, int n, int d){
  NumericMatrix Mult_h(n,m);
  // NumericMatrix Mult_C(n,m);
  NumericVector norm_t(n);

  // NumericMatrix U_hC(m,2*n-2);
  NumericMatrix U_hC(m,n-1);

  for (int t=0; t<m; t++){

    for (int i=0; i<n; i++){
      norm_t(i)=norm(t*n + i);
    }

    NumericVector A_norm_h = vecmult(A_h, norm_t);
    // NumericVector A_norm_C = vecmult(A_C, norm_t);
    for (int i=0; i<A_norm_h.length(); i++){
      Mult_h(i,t)=A_norm_h(i);
      // Mult_C(i,t)=A_norm_C(i);
    }

    NumericMatrix h_Mult(d*(n-1), n);
    // NumericMatrix C_Mult(d*(n-1),n);
    NumericMatrix h_Mult_sum(d, n-1);
    // NumericMatrix C_Mult_sum(d, n-1);

    parallel_for(d, [&](int start, int end){
      for(int D=start; D<end; ++D){
        for(int i=0; i<n-1; i++){
          for(int j=i+1; j<n; j++){  // add multiplier
            h_Mult(D*(n-1)+i,j)=hC_Obs(D*(n-1)+i,j)*(Mult_h(i,t)+Mult_h(j,t));  //h
            // C_Mult(D*(n-1)+i,j)=hC_Obs(D*(n-1)+i,j+n)*(Mult_C(i,t)+Mult_C(j,t));  //C
          }
        }

        for(int i=0;i<n-1; i++){
          for(int j=0;j<=i; j++){
            for(int k=i+1;k<n; k++){
              h_Mult_sum(D,i) += h_Mult(D*(n-1)+j,k);
              // C_Mult_sum(D,i) += C_Mult(D*(n-1)+j,k);
            }
          }
        }

      }
    });

    for (int k=0; k<n-1; k++){
      double help_h =0;
      // double help_C =0;
      for(int D=0; D<d; D++){
        help_h += h_Mult_sum(D,k)*h_Mult_sum(D,k);
        // help_C += C_Mult_sum(D,k)*C_Mult_sum(D,k);
      }
      // U_{n,k}^(t) for Wilcoxon and Spatial Sign, saved in one matrix.
      //    Each row corresponds to one bootstrap iteration
      U_hC(t,k)=std::sqrt(help_h);  //h
      // U_hC(t,k+(n-1))=std::sqrt(help_C);  //C
    }
  }

  return U_hC;
}


// matrix-vector multiplication
// [[Rcpp::export]]
NumericVector vecmult(NumericMatrix M, NumericVector v){
  int r = M.nrow();
  int c = M.ncol();

  NumericVector w(r);
  for(int i=0; i< r;i++) {
    for(int k=0;k<c;k++){
      w(i) += M(i,k) * v(k);
    }

  }
  return w;
}


// [[Rcpp::export]]
NumericMatrix outerProd(NumericVector v, NumericVector w){
  int n=v.length();
  int m=w.length();
  NumericMatrix M(n,m);
  for(int i=0;i<n;i++){
    M(i,_) = v(i) * w;
  }
  return M;
}


// [[Rcpp::export]]
NumericVector find_rowmax(NumericMatrix U_hC){
  NumericVector max_hC(U_hC.nrow());//2);
  // NumericVector vC(n-1);
  for (int t=0; t<U_hC.nrow(); t++){
    max_hC(t) = max(U_hC(t,_));
    // max_hC(t,1) = max(vC);
  }
  return max_hC;
}


// Calculate U_{n,k} for k=1, .., n-1
// [[Rcpp::export]]
NumericMatrix fill_T(NumericMatrix hC_Obs, int n, int d){
  NumericMatrix h_Obs_sum(d, n-1);
  // NumericMatrix C_Obs_sum(d, n-1);

  for(int D=0; D<d; D++){
    for(int i=0;i<n-1; i++){
      for(int j=0;j<=i; j++){
        for(int k=i+1;k<n; k++){
          h_Obs_sum(D,i) += hC_Obs(D*(n-1)+j,k);
          // C_Obs_sum(D,i) += hC_Obs(D*(n-1)+j,k+n);
        }
      }
    }
  }

  NumericMatrix T_hC((n-1),1);//2);
  for (int k=0; k<n-1; k++){
    double help_h =0;
    // double help_C =0;
    for(int D=0; D<d; D++){
      help_h += h_Obs_sum(D,k) * h_Obs_sum(D,k);
      // help_C += C_Obs_sum(D,k)*C_Obs_sum(D,k);
    }
    T_hC(k,0)=std::sqrt(help_h);
    // T_hC(k,1)=std::sqrt(help_C);
  }

  return T_hC;
}


// Calculate h(v1,v2) = S(v1-v2), S( ) spatial sign function
// [[Rcpp::export]]
NumericVector h_cpp(NumericVector v1, NumericVector v2){
  int n=v1.length();
  NumericVector res(n);

  double help=0;
  int eq=0;
  for (int i=0; i<n; i++){
    if (v1(i)==v2(i)) eq += 1;
    help += (v1(i)-v2(i)) * (v1(i)-v2(i));
  }

  if( eq != n) {
    for (int i=0; i<n; i++){
      res(i)=(v1(i)-v2(i)) / std::sqrt(help);
    }
  }

  return res;
}
