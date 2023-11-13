
compute_Tn1 <- function(X, M = 100000, W = NULL, space = "BM", ...) {
  n <- ncol(X)

  if (is.null(W)) {
    W <- computeSpaceMeasuringVectors(M = M, X = X, space = space)
  } else {
    M <- ncol(W)
  }

  Rcpp::cppFunction('double dot_integrate1(NumericVector v) {
    return sum(v[Rcpp::Range(1, v.size()-1)] + v[Rcpp::Range(0, v.size()-2)]) / (2 * (v.size() - 1));
  }')

  intVal <- sapply(1:M, function(v, W, X1, n) {
    dot_integrate1(abs(.Zn(W[, v], X1))^2)
  }, W = W, X1 = X, n = n)

  1 / M * sum(intVal)
}

compute_Tn2 <- function(X, M = 100000, W = NULL, space = "BM", ...) {
  n <- ncol(X)

  if (is.null(W)) {
    W <- computeSpaceMeasuringVectors(M = M, X = X, space = space)
  } else {
    M <- ncol(W)
  }

  intVal <- sapply(1:M, function(v, W, X1, n) {
    dot_integrate(abs(.Zn(W[, v], X1))^2)
  }, W = W, X1 = X, n = n)

  1 / M * sum(intVal)
}
