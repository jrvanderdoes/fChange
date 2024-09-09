#' Characteristic Functions Change Point Detection - Welch Approximation
#'
#' @param X funts object
#' @param M Numeric for the number of vectors to measure space
#' @param J Numeric for the level of discretion
#' @param h Numeric for bandwidth
#' @param K Function for bandwidth
#' @param space String for the space to measure data
#' @param alpha Numeric for significance
#'
#' @return List with the following items:
#'    'cutoff': Cutoff calculated for Tn using Welch approximation,
#'    'value': Test statistic Tn for the data
#'    'detected': Boolean if Tn is significant or not
#'    'alpha': Significance of the result
#' @export
#'
#' @examples
#' characteristic_change_welch(funts(electricity))
characteristic_change_welch <- function(X,
                                     M = 20, J=50,
                                     h = ncol(X$data)^(1 / 3),
                                     K = bartlett_kernel,
                                     space = "BM",
                                     alpha = 0.05){
  X <- .check_data(X)

  # Generate Noise
  W <- computeSpaceMeasuringVectors(M = M, X = X, space = space)

  # Determine Number of Iterations
  val_cutoff <- .compute_Welch(X, alpha = alpha, W = W,
                         M = M, J=J, h = h, K = K )
  val_Tn <- compute_Tn(X=X$data, W=W, J=J)

  list(
    "cutoff" = val_cutoff,
    "value" = val_Tn,
    'detected' = val_Tn > val_cutoff,
    'alpha'=alpha
  )
}


#' Welch Approximation to Tn
#'
#' This function approximates the Tn statistic using the Welch approximation.
#'
#' @param X Numeric data.frame with evaled points on rows and fd objects in columns
#' @param alpha (Optional) Numeric value in (0, 1) for significance. Default is 0.05
#' @param TVal (Optional) Numeric value indicate the number of FDs. Default is the
#'     number of columns in X.
#' @param W (Optional) Data.frame for the functions to integrate against
#'     used for both muHat and sigmaHat2. Default of NULL will use Brownian
#'     motion.
#' @param W1 (Optional)  data.frame for the functions to integrate against
#'     used for sigmaHat2. Default of NULL will use Brownian motion.
#' @param M (Optional) Numeric value for number of functions to integrate against
#'     (ncol(W)). Default is 100.
#' @param h (Optional) Numeric value indicating the amount of blocking. Default
#'     is Tval^(1/3)
#' @param K (Optional) Function used for CHat, indicating Kernel. Default is the
#'     Bartlett Kernel.
#' @param ... Unused parameters, but added for use in other functions.
#'
#' @return Numeric value indicating cutoff from Welsh approximation
#'
#' @examples
#' .compute_Welch(X=electricity[,1:40], h = 0)
.compute_Welch <- function(X, alpha = 0.05,
                     W = NULL, M = 100, J = 50,
                     h = ncol(funts(X)$data)^(1 / 3),
                     K = bartlett_kernel ) {
  x <- .check_data(X)

  if (is.null(W)) {
    W <- computeSpaceMeasuringVectors(M, "BM", X)
  }

  m_gamma <- 0
  for(i in 1:M){
    Cre <- .estimD(K=K, h=h, X=X$data,
                   lfun=cos, v=W[,i], lfunp=cos, vp=W[,i])
    Cim <- .estimD(K=K, h=h, X=X$data,
                   lfun=sin, v=W[,i], lfunp=sin, vp=W[,i])

    m_gamma <- m_gamma + ( Cre + Cim )
  }
  m_gamma <- m_gamma / (6*M)

  sigma2_gamma <- sigma2_gamma1 <- 0
  for (i in 1:M) {
    val <- val1 <- 0
    for (j in 1:M){
      if (i==j) next
      Cre <- .estimD(K=K, h=h, X=X$data,
                     lfun=cos, v=W[,i], lfunp=cos, vp=W[,j])
      Cri <- .estimD(K=K, h=h, X=X$data,
                     lfun=cos, v=W[,i], lfunp=sin, vp=W[,j])
      Cim <- .estimD(K=K, h=h, X=X$data,
                     lfun=sin, v=W[,i], lfunp=sin, vp=W[,j])

      val <- val + (Cre^2+2*Cri^2+Cim^2)
    }
    val <- val / (M-1)
    sigma2_gamma <- sigma2_gamma + val
  }
  sigma2_gamma <- (2 / 90) * sigma2_gamma / M

  beta <- sigma2_gamma / (2*m_gamma)
  nu <- (2*m_gamma^2) / sigma2_gamma

  beta * stats::qchisq(1 - alpha, df = nu)
}
