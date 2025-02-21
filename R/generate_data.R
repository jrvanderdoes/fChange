#' Generate Functional Data
#'
#' This function allows generation of functional data according to several
#'  approaches: \code{bbridge}, \code{bmotion}, \code{kl}, \code{ou}, and \code{far1}.
#'
#' @param data_details List of named lists indicating parameters for each data group.
#'  Each process can use different parameters, given below.
#'  \itemize{
#'    \item **bmotion**: Brownian motion contains:
#'      \itemize{
#'        \item **N**: Numeric indicating the number of observations.
#'        \item **sd**: Numeric for the standard deviation of the observations.
#'      }
#'    \item **bbridge**: Brownian bridge contains:
#'      \itemize{
#'        \item **N**: Numeric indicating the number of observations.
#'        \item **sd**: Numeric for the standard deviation of the observations.
#'      }
#'    \item **kl**: Karhunen-Loeve expansion contains:
#'      \itemize{
#'        \item **N**: Numeric indicating the number of observations.
#'        \item **distribution**: Distribution of the errors. Options include:
#'          \code{binomial}, \code{exponential}, \code{laplace}, \code{normal},
#'          \code{t} (add dof argument), \code{gamma} (add shape argument),
#'           and \code{cauchy}.
#'        \item **eigenvalues**: Numerics for the eigenvalues of the given
#'          distribution (value for each in basis).
#'        \item **mean**: Numeric for the mean of the group.
#'        \item **dependence**: Strength of dependence between observation.
#'        \item **basis**: fda \code{basisfd} object.
#'        \item **dof**: (Optional) Numeric for the degrees of freedom if using
#'          a \code{t} distribution.
#'        \item **shape**: (Optional) Numeric for the shape if using
#'          a \code{gamma} distribution.
#'      }
#'    \item **ou**: Ornstein-Uhlenbeck process requires:
#'      \itemize{
#'        \item **N**: Numeric indicating the number of observations.
#'        \item **dependence**: Strength of dependence between observation.
#'      }
#'    \item **far1**: Functional autoregressive process of order 1 contains:
#'      \itemize{
#'        \item **N**: Numeric indicating the number of observations.
#'        \item **dependence**: Strength of dependence between observation.
#'        \item **sd**: Numeric for the standard deviation of the observations.
#'        \item **vary**: Boolean if the starting value each observation should
#'          be 0 or vary. It does this by dropping first value. If resolution is
#'          given as a number, it can adjust so that the length is the same. If
#'          resolution is a vector, the resolution will be one smaller.
#'      }
#'  }
#'
#' @return dfts object for the generated data
#' @export
#'
#' @examples
#' result <- generate_data(
#'   resolution=15,
#'   data_details =list('bmotion'=list('N'=100, 'sd'=1),
#'                      'bbridge'=list('N'=100, 'sd'=1),
#'                      'bbridge'=list('N'=100, 'sd'=1),
#'                      'kl'=list('N'=100,
#'                                'distribution'='Normal',
#'                                'eigenvalues'=1/1:4,
#'                                'mean'=0, 'dependence'=0,
#'                                'basis'=fda::create.bspline.basis(),
#'                                'sd'=1),
#'                      'ou'=list('N'=100, 'dependence'=0 ) ,
#'                      'far1'=list('N'=100, 'dependence'=0,
#'                                  'sd'=1,'vary'=FALSE) )
#' )
generate_data <- function(resolution, data_details, burnin=100){
  # Setup and Check
  resolution_base = resolution
  if(length(resolution)==1){
    resolution <- seq(0,1, length.out=resolution)
  }

  poss_process <- c('bbridge','bmotion','kl','ou','far1')
  data <- matrix()
  list_names <- names(data_details)

  # Go through each defined groups
  for(i in 1:length(list_names)){
    # Parameters
    process <- .verify_input(list_names[i], poss_process)
    N_data <- data_details[[i]]$N
    if(i==1){
      N <- N_data + burnin
    } else{
      N <- N_data
    }

    # Generate Data
    data_tmp <- switch(process,
                       bbridge = {
                         results_tmp <- generate_brownian_bridge(N = N, v = resolution,
                                                                 sd = data_details[[i]]$sd)
                         results_tmp$data[, (N-N_data+1):N]
                       },
                       bmotion = {
                         results_tmp <- generate_brownian_motion(N = N, v = resolution,
                                                                 sd = data_details[[i]]$sd)
                         results_tmp$data[, (N-N_data+1):N]
                       },
                       kl = {
                         results <- generate_karhunen_loeve(N = N,
                                                            eigenvalues = data_details[[i]]$eigenvalues,
                                                            basis = data_details[[i]]$basis,
                                                            means = data_details[[i]]$mean,
                                                            distribution = data_details[[i]]$distribution,
                                                            resolution = resolution,
                                                            dependence = data_details[[i]]$dependence,
                                                            burnin = N-N_data,
                                                            silent = TRUE,
                                                            prev_eps = data_details[[i]]$prev_eps,
                                                            dof=data_details[[i]]$dof,
                                                            shape=data_details[[i]]$shape)
                         if(i < length(data_details))
                           data_details[[i+1]]$prev_eps <- results$prev_eps

                         results$data[, (N-N_data+1):N]
                       },
                       ou = {
                         results <- generate_ornstein_uhlenbeck(resolution = resolution, N = N,
                                                 rho = data_details[[i]]$dependence)

                         results$data[, (N-N_data+1):N]
                       },
                       far1 = {
                         if(length(resolution_base)==1 & data_details[[i]]$vary){
                           resolution1 <- seq(0, 1, length.out=resolution_base+1)
                         }else if(data_details[[i]]$vary){
                           stop(paste0('A specified resolution is given;',
                                       'however, with paramter `vary=TRUE`, we must lose',
                                       'a resolution point. It is unclear how to do so.',
                                       'Specify resolution as a number or do not vary'))
                         } else{
                           resolution1 <- resolution
                         }

                         results <- generate_far1(N = N, resolution = resolution1,
                                                  sd=data_details[[i]]$sd,
                                                  dependence = data_details[[i]]$dependence,
                                                  dropFirst=data_details[[i]]$vary)

                         results$data[, (N-N_data+1):N]
                       }
    )

    if(i==1){
      data <- data_tmp
    }else{
      data <- cbind(data,data_tmp)
    }
  }

  dfts(data, intratime = resolution, labels = 1:ncol(data))
}
