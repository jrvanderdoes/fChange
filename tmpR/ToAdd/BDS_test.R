#' Nonlinear independence BDS test
#'
#' @param X funts object
#' @param m embedding dimension hyperparameter. Common choices are 2, 3 or 4.
#'  Default is 3.
#' @param r distance threshold value. Common choices of r are 0.75s.d.(X),
#'  s.d.(X) and 1.25s.d.(X). Default is sd(X$data)
#'
#' @return List with two elements
#' \itemize{
#'  \item statistic: BDS statistic for the data
#'  \item pvalue: Computed p-value
#' }
#' @export
#'
#' @references Huang, Xin, Han Lin Shang, and Tak Kuen Siu. “A Nonlinearity
#'  and Model Specification Test for Functional Time Series,” 2023.
#'  https://doi.org/10.48550/arxiv.2304.01558.
#'
#' @examples
#' BDS_test(funts(electricity))
#' BDS_test(generate_brownian_bridge(500))
#' BDS_test(generate_brownian_motion(500))
#' BDS_test(funts(generateFAR1(20,500)))
#' BDS_test(funts(generateFAR1(20,500,d=0)))
#' BDS_test(generateOU(20,500))
BDS_test <- function(X, m=3, r=sd(X$data)){
  n <- ncol(X$data)

  indicator <- matrix(NA, nrow=n, ncol=n)
  for (i in 1:n){
    for (j in i:n){
      indicator[j,i] <- indicator[i,j] <-
        sqrt(sum( (X$data[,i]-X$data[,j])^2 ) ) < r
    }
  }
  bds_c_T <- ( sum(rowSums(indicator)^2)-3*sum(indicator)+2*n ) /
    ( n*(n-1)*(n-2) )

  counts <- (n-m+1)*(n-m) / 2
  M <- indicator[1:(n-m+1), 1:(n-m+1)]
  for (m_ind in 2:m){
    M <- M * indicator[m_ind:(n+m_ind-m), m_ind:(n+m_ind-m)]
  }
  bds_c_m <- sum( M[upper.tri(M, diag = FALSE)] ) / counts

  M <- indicator[m:n,m:n]
  bds_c_m1 <- sum( M[upper.tri(M, diag = FALSE)] ) / counts

  counts <- n*(n-1)/2
  M <- indicator[1:n,1:n]
  C <- sum( M[upper.tri(M, diag = FALSE)] ) / counts

  A <- 0
  for(j in 1:(m-1)){
    A <- A + ( bds_c_T^(m-j) * C^(2*j) )
  }
  sigma <- 4 * (bds_c_T^m + 2*A + (m-1)^2*C^(2*m) - m^2*bds_c_T*C^((2*m)-2))
  BDS_statistic <- ( bds_c_m-bds_c_m1^m ) * sqrt(n-m+1) / sqrt(sigma)
  p_value <- min(c(pnorm(BDS_statistic), 1-pnorm(BDS_statistic)))

  list('statistic'=BDS_statistic, 'pvalue'=p_value)
}
