#' Generate functional data
#'
#' \code{generate_karhunen_loeve} generates functional data via KL expansion.
#' This can include change points in any combination of the following:
#' \itemize{
#'   \item Mean
#'   \item Distribution
#'   \item Eigenvalue(s)
#'   \item Eigenvector(s)
#' }
#' In this sense, the function creates m 'groups' of discretely observed
#'  functions with similar properties.
#'
#' @param N Vector of Numerics. Each value in N is the number of observations
#'  for a given group.
#' @param eigenvalues Vector of eigenvalues. Length 1 or m.
#' @param basis A list of bases (eigenfunctions), length m.
#' @param means A vector of means, length 1 or N.
#' @param distribution A vector of distributions, length 1 or m.
#' @param intratime A vector of points indicating the points to evaluate the
#'     functions on.
#' @param dependence Numeric \[0,1\] indicating strength of VAR(1) process.
#' @param burnin A numeric value indicating the number of burnin trials.
#' @param silent A Boolean that toggles running output.
#' @param shape Numeric for degrees of freedom with t-distribution
#' @param dof Numeric for shape with gamma distribution (rate is set to 1)
#' @param prev_eps previous epsilon for dependence across groups.
#'
#' @return List with (1) data (N-by-m) and (2) previous errors.
#' @export
#'
#' @examples
#' dat1 <- generate_karhunen_loeve(
#'   N=100, eigenvalues=c(1/(1:3)), basis=fda::create.bspline.basis(nbasis=3,norder=3),
#'   means=0, distribution='Normal',
#'   intratime=seq(0,1,0.1), dependence=0, burnin=100, silent=TRUE, dof=NULL, shape=NULL,
#'   prev_eps=NULL)
#' dat2 <- generate_karhunen_loeve(
#'   N=50, eigenvalues=c(1/(1:4)), basis=fda::create.bspline.basis(nbasis=4),
#'   means=5, distribution='exponential',
#'   intratime=seq(0,1,0.1), dependence=0, burnin=100, silent=TRUE, dof=NULL, shape=NULL,
#'   prev_eps=dat1$prev_eps)
#'
#' dat <- dfts(cbind(dat1$data$data, dat2$data$data),intratime = dat1$data$intratime)
generate_karhunen_loeve <- function(
    N, eigenvalues, basis, means, distribution,
    intratime, dependence=0, burnin=100, silent=TRUE, dof=NULL, shape=NULL,
    prev_eps=NULL) {
  ## Verification
  m <- length(basis$names)

  .check_length <- function(n, vecs, name){
    if(length(vecs)==1){
      vecs <- rep(vecs, n)
    }else if(length(vecs)!= n){
      stop(paste0('Not enough ',name,' given'))
    }

    vecs
  }

  eigenvalues <- .check_length(m, eigenvalues, 'eigenvalues')
  distribution <- .check_length(m, distribution, 'distribution')

  means <- .check_length(N, means, 'means')

  ## Prepare to generate

  # Setup psi
  .getPsi <- function(m, eigenvalues, dependence) {
    psi <- list()
    normsSD <- stats::rnorm(m, mean = 0, sd = 1)

    groupSD <- normsSD[1:m] * sqrt(eigenvalues)
    psi0 <- groupSD %*% t(groupSD)
    psi0 <- psi0 / sqrt(sum(psi0^2)) ## TODO:: Check this

    dependence * psi0
  }

  psi <- .getPsi(m, eigenvalues, dependence)

  # Burnin for VAR (if not given)
  if(is.null(prev_eps)){
    prev_eps <- data.frame(matrix(0, ncol = length(intratime), nrow = m))
    if(burnin>0){
      for (j in 1:burnin) {
        waste <- .KL_Expansion(
          eigenvalues = eigenvalues,
          basis = basis,
          means = means[1],
          distribution = distribution,
          resolution = intratime,
          prev_eps = prev_eps,
          psi = psi,
          dof = dof,
          shape = shape
        )

        prev_eps <- waste[[2]]
      }
    }
  }

  # If Num of Eigs increases or decreases (only at changes)
  psiDim1 <- dim(psi)[1]
  pepDim1 <- dim(prev_eps)[1]

  if (psiDim1 != pepDim1) {
    if (psiDim1 > pepDim1) {
      # Bind row of 0s to the bottom if didn't have a value previously
      for (k in 1:(psiDim1 - pepDim1)) {
        prev_eps <- rbind(prev_eps, 0)
      }
    } else if (psiDim1 < pepDim1) {
      # Drop Rows if not needed
      for (k in 1:(pepDim1 - psiDim1)) {
        prev_eps <- prev_eps[-nrow(prev_eps), ]
      }
    }
  }

  # Generate Data
  data <- data.frame(matrix(NA, ncol = N, nrow = length(intratime)))
  for (j in 1:N) {
    result <- .KL_Expansion(
      eigenvalues = eigenvalues,
      basis = basis,
      means = means[j],
      distribution = distribution,
      resolution = intratime,
      prev_eps = prev_eps,
      psi = psi,
      dof = dof,
      shape = shape
    )

    data[, j] <- result[[1]]
    prev_eps <- result[[2]]
  }

  list('data'= dfts(data, intratime = intratime, labels = 1:ncol(data)),
       'prev_eps'=prev_eps)
}


#' KL Expansion Computation
#'
#' @inheritParams generate_karhunen_loeve
#' @param psi Matrix for dependence with previous errors.
#'
#' @return List with (1) data (N-by-m) and (2) previous errors.
#'
#' @noRd
#' @keywords internal
.KL_Expansion <- function(
    eigenvalues, basis, means, distribution, resolution, prev_eps, psi, shape, dof, ...)
  {
  # Setup
  n <- length(resolution)
  D <- length(eigenvalues)

  # Generate - No Loop
  eval_basis <- fda::eval.basis(resolution, basis)

  xi <- sapply(1:length(eigenvalues), function(e, eigenvalues, distribution, dof, shape, ...) {
    .KL_errors(distribution = distribution[e], sd = sqrt(eigenvalues[e]),
               n = 1, shape=shape, dof=dof, ...)
    },
  eigenvalues = eigenvalues, distribution = distribution, dof = dof, shape=shape, ...
  )

  Zeta <- tryCatch(eval_basis * xi[col(eval_basis)],
                   error = function(e) {
                     stop(call. = F, paste0(
                       "Check number of eigenvalues given. ",
                       "It does not match number of basis functions. ",
                       "Note, did you account for the constant function if ",
                       "it is in the basis?"
                     ))
                   }
  )

  eps <- Zeta + t(psi %*% as.matrix(prev_eps))
  X <- means + rowSums(eps)

  list(X, t(eps))
}


#' KL Expansion Errors
#'
#' @inheritParams generate_karhunen_loeve
#' @param sd Standard deviation of errors (sqrt(eigenvalues))
#' @param n Numeric number of errors to generate
#'
#' @return Errors as a vector
#'
#' @noRd
#' @keywords internal
.KL_errors <- function(distribution, sd, n = 1, dof = NULL, shape = 1, rate = 1) {
  ## This function give centered distributions with eig^2 var

  if (tolower(distribution) == "normal") {
    xi <- stats::rnorm(n, mean = 0, sd = sd)
  } else if (tolower(distribution) == "binomial") {
    if (sd == 0) {
      return(rep(0, n))
    }

    mean_val <- 10 * sd^2 # arbitrary, must exceed var
    p <- 1 - sd^2 / mean_val
    size <- round(mean_val / p)

    xi <- stats::rbinom(n = n, size = size, p = p) - mean_val
  } else if (tolower(distribution) == "exponential") {
    xi <- stats::rexp(n, rate = 1 / sd) - sd
  } else if (tolower(distribution) == "t") {
    xi <- stats::rt(n, dof)
    if(dof>2)
      xi <- xi * sqrt(sd^2 * (dof - 2) / dof)
  } else if (tolower(distribution) == "cauchy") {
    xi <- stats::rt(n, 1)
  }  else if (tolower(distribution) == "gamma") {
    xi <- ( stats::rgamma(n,shape = shape, rate = rate) - shape/rate ) /
      sqrt( shape * (1/rate)^2 ) * sd
  } else if (tolower(distribution) == "laplace") {
    if (!requireNamespace("jmuOutlier", quietly = TRUE)) {
      stop(paste0("Please install 'jmuOutlier'."), call. = FALSE)
    }
    xi <- jmuOutlier::rlaplace(n, mean = 0, sd = sd)
  } else {
    stop(paste("Sorry, distribution", distribution, "not implemented yet"))
  }

  xi
}
