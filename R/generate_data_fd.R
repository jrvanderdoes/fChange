#' Generate functional data
#'
#' \code{generate_data_fd} generates functional data via KL expansion.
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
#' @param ns A numerical vector of length m
#'
#'     Indicates the number of functional objects (groups) created using first
#'     set of parameters
#'
#' @param eigsList A list of vectors of length $1$ or $m$ with the eigenvalues
#'     for each group
#' @param basesList A list of bases (eigenfunctions), length 1 or m
#'
#'     Define the basis using fda on c(0, 1) to ensure it works
#'     (TODO:: Remove this restriction)
#' @param meansList A list of means, length 1 or m, for each group
#' @param distsArray A vector of distributions, length 1 or m, for each group
#' @param evals A vector of points indicating the points to evaluate the
#'     functions on
#' @param kappasArray A vector of Kappas, length 1 or m, for the strength of the
#'    VAR(1) process
#' @param burnin A numeric value indicating the number of burnin trials
#'    This is only necessary when kappa>0
#' @param silent A Boolean that toggles running output
#' @inheritDotParams .generateX dof
#'
#' @return A data.frame of m columns length(evals) rows (TODO:: Verify)
#'
#' @export
#'
#' @examples
#' # Create 200 functions with a midway change point. The change point
#' #     is a change point in the eigenvalues, eigenfunctions, means,
#' #     distributions, and VAR(1) strength
#' data_KL <- generate_data_fd(
#'   ns = c(25, 25),
#'   eigsList = list(
#'     c(3, 2, 1, 0.5),
#'     c(2, 3, 2)
#'   ),
#'   basesList = list(
#'     fda::create.bspline.basis(nbasis = 4, norder = 4),
#'     fda::create.fourier.basis(nbasis = 2)
#'   ),
#'   meansList = c(0, 0.5),
#'   distsArray = c("Normal", "Binomial"),
#'   evals = seq(0, 1, 0.05),
#'   kappasArray = c(0, 0.5)
#' )
generate_data_fd <- function(ns,
                             eigsList,
                             basesList,
                             meansList,
                             distsArray,
                             evals,
                             kappasArray = c(0),
                             burnin = 100,
                             silent = FALSE) {
  # ns is a vector with length m for the number of data runs until next CP
  # - i.e. c(10,10,10) has 10 length TS then CP followed by 10 and another CP
  # eigsList is a list of vectors giving the eigenvalues for each distribution
  # - Length of list be 1 or m
  # basesList is a list of bases (eigf) for each distribution
  # - Length should be 1 or m
  # - Define basis on (0, 1)
  # meansList is a list of means for each distribution
  # - Length should be 1 or m
  # distsArray is a vector of distributions to each run
  # - Length should be 1 or m
  # evals is the vector of points to evaluate
  # kappasArray is a vector of Kappas for VAR
  # - Length should be 1 or m

  ## Verification
  m <- length(ns)

  eigsList <- .checkLength(eigsList, "eigsList", m)
  basesList <- .checkLength(basesList, "basesList", m)
  meansList <- .checkLength(meansList, "meansList", m)
  distsArray <- unlist(.checkLength(distsArray, "distsArray", m))
  kappasArray <- unlist(.checkLength(kappasArray, "kappaArray", m))

  # Prepare to generate
  data <- data.frame(matrix(NA, ncol = sum(ns), nrow = length(evals)))
  addIdx <- 0

  # Setup psi
  Ds <- 1:m
  for (i in 1:m) {
    Ds[i] <- length(eigsList[[i]])
  }
  psi <- .getPsiList(Ds, eigsList, kappasArray)

  # Burnin for VAR
  peps <- data.frame(matrix(0, ncol = length(evals), nrow = Ds[1]))
  for (j in 1:burnin) {
    waste <- .generateData_KL_Expansion(
      eigs = eigsList[[1]],
      basis = basesList[[1]],
      means = meansList[[1]],
      dist = distsArray[1],
      evals = evals,
      peps = peps,
      psi = psi[[1]],
      ...
    )

    peps <- waste[[2]]
  }

  for (i in 1:m) {
    if (!silent) cat(paste0("Running setup ", i, "/", m, "\n"))

    # If Num of Eigs increases or decreases (only at CPs)
    psiDim1 <- dim(psi[[i]])[1]
    pepDim1 <- dim(peps)[1]
    if (psiDim1 != pepDim1) {
      if (psiDim1 > pepDim1) {
        # Bind row of 0s to the bottom if didn't have a value previously
        for (k in 1:(psiDim1 - pepDim1)) {
          peps <- rbind(peps, 0)
        }
      } else if (psiDim1 < pepDim1) {
        # Drop Rows if not needed
        for (k in 1:(pepDim1 - psiDim1)) {
          peps <- peps[-nrow(peps), ]
        }
      }
    }

    for (j in 1:ns[i]) {
      result <- .generateData_KL_Expansion(
        eigs = eigsList[[i]],
        basis = basesList[[i]],
        means = meansList[[i]],
        dist = distsArray[i],
        evals = evals,
        peps = peps,
        psi = psi[[i]],
        ...
      )

      data[, addIdx + j] <- result[[1]]
      peps <- result[[2]]
    }
    addIdx <- addIdx + ns[i] # j
  }

  data
}


#' Check Length
#'
#' This (internal) function checks the length of inputs to verify all is well
#'  before spending computational time.
#'
#' @param dataList List (or vector at times) containing the data.
#' @param name String of variable being checked so error can be informative.
#' @param m Numeric integer indicating the number of segments that will be
#'  built, i.e. the data should generally be this length or 1 (if all the same).
#'
#' @return dataList, perhaps extending it as needed (making repeats).
#'
#' @noRd
.checkLength <- function(dataList, name, m) {
  if (length(dataList) == 1) {
    if (fda::is.basis(dataList[[1]])) {
      retList <- list()
      for (i in 1:m) {
        retList <- append(retList, dataList)
      }
      dataList <- retList
    } else {
      dataList <- list(rep(dataList[[1]], m))
    }
  } else if (length(dataList) != m) {
    stop(paste(name, "is length", length(dataList), "not 1 or", m, "\n"))
  }

  dataList
}


#' Setup psi
#'
#' This (internal) function sets up the psi for dependence of each segment.
#'
#' @param D Vector of numerics for the number of eigenvalues in each segment
#' @inheritParams generate_data_fd
#'
#' @return List with the dependence for each segment
#'
#' @noRd
.getPsiList <- function(D, eigsList, kappasArray) {
  psi <- list()
  normsSD <- stats::rnorm(max(D), mean = 0, sd = 1)

  for (i in 1:length(D)) {
    groupSD <- normsSD[1:D[i]] * sqrt(eigsList[[i]])
    psi0 <- groupSD %*% t(groupSD)
    psi0 <- psi0 / sqrt(sum(psi0^2)) ## TODO:: Check this
    psi[[i]] <- kappasArray[i] * psi0
  }

  psi
}


#' Generate Data - KL Expansion
#'
#' This (internal) function performs the KL expansion to create FD object.
#'
#' @param eigs Vector of numeric values indicating the eigenvalues of interest
#' @param basis FDA basis object
#' @param means Numeric for mean of FD object
#' @param dist String indicating the distribution to use. Options include Normal,
#'             Binomial, Exponential, and t.
#' @param evals A vector of points indicating the points to evaluate the
#'     functions on
#' @param peps Vector of numerics for the previous epsilons
#' @param psi Numeric for psi value.
#'
#' @return List of X, observed points in FD, and eps, the epsilons for next
#'         observation.
#'
#' @noRd
.generateData_KL_Expansion <- function(eigs, basis, means, dist,
                                       evals, peps, psi, ...) {
  # Setup
  n <- length(evals)
  D <- length(eigs)

  # X <- rep(0, n)
  # Zeta <- eps <-
  #   data.frame(matrix(NA,ncol=n,nrow=D)) # Matrix with col as time, row as dimension

  # Verify
  if (length(means) == 1) {
    means <- rep(means, n)
  } else if (length(means) != n) {
    stop(paste("Length of means is", length(means), "not 1 or", n))
  }

  # Generate - No Loop
  eval_basis <- fda::eval.basis(evals, basis)
  # Row for each time, columns for eigen
  xi <- sapply(eigs, function(e, dist, n) {
    .generateXi(dist = dist, sd = sqrt(e), n = n, ..)
  },
  dist = dist, n = n
  )

  Zeta <- tryCatch(xi * eval_basis,
    error = function(e) {
      stop(call. = F, paste0(
        "Check number of eigenvalues given. ",
        "It does not match number of basis functions."
      ))
    }
  )

  eps <- Zeta + t(psi %*% as.matrix(peps))
  X <- means + rowSums(eps)

  # Generate
  # for(t in 1:n){
  #   eval_basis <- fda::eval.basis(evals[t], basis)
  #   for(j in 1:D){
  #     xi <- .generateXi(dist=dist, sd=sqrt(eigs[j]))
  #     xi2[t,j] <- xi
  #     Zeta[j,t] <- xi * eval_basis[j]
  #   }
  #
  #   eps[,t] <- Zeta[,t] + psi %*% peps[,t]
  #   X[t] <- means[t] + sum(eps[,t])
  # }

  list(X, t(eps))
}


#' Generate Xi
#'
#' This (internal) function generates the noise for KL expansion terms
#'
#' @param dist String for distribution of interest. Options are: Binomial,
#'  Exponential, laplace, Normal, and t. Upcoming is cauchy.
#' @param sd Numeric for the standard deviation of the variable
#' @param n (Optional) Numeric for the number of observations. Default is 1.
#' @param dof (Optional) Numeric value indicating the degrees of freedom used
#'  in the t-distribution. Default is 3. Note, if value is less than 3, the sd
#'  parameter cannot be followed.
#'
#' @return Vector of numerics for observations of the variable.
#'
#' @noRd
.generateXi <- function(dist, sd, n = 1, dof=NULL) {
  ## This function give centered distributions with eig^2 var

  if (dist == "Normal") {
    xi <- stats::rnorm(n, mean = 0, sd = sd)
  } else if (dist == "Binomial") {
    if (sd == 0) {
      return(rep(0, n))
    }

    mean <- 10 * sd^2 # arbitrary, must exceed var
    p <- 1 - sd^2 / mean
    size <- round(mean / p)

    xi <- stats::rbinom(n = n, size = size, p = p) - mean
  } else if (dist == "Exponential") {
    xi <- stats::rexp(n, rate = 1 / sd) - sd
  } else if (dist == "t") {
    xi <- stats::rt(n, dof)
    if(dof>2)
      xi <- xi * sqrt(sd^2 * (dof - 2) / dof)
  } else if (dist == "cauchy") {
    stop("Sorry Problem with cauchy directly. Try t-distribution")
    xi <- stats::rcauchy(n)
  } else if (dist == "laplace") {
    if (!requireNamespace("jmuOutlier", quietly = TRUE)) {
      stop(paste0("Please install 'jmuOutlier'."))
    }
    xi <- jmuOutlier::rlaplace(n, mean = 0, sd = sd)
  } else {
    stop(paste("Sorry, dist", dist, "not implemented yet"))
  }

  xi
}
