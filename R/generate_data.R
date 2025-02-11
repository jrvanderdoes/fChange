#' Generate Functional Data
#'
#' This function allows generation of functional data according to several
#'  approaches: \code{bbridge}, \code{bmotion}, \code{kl}, \code{ou}, and \code{far1}.
#'
#' @param general List of general parameters for entire data sequence.
#'  \itemize{
#'    \item **resolution**: Numeric(s). Either a number indicating the
#'      observational resolution which are evenly spaced on \[0,1\] or a vector of
#'      points which indicate the evalution points. Resolution does not change
#'      between generated data.
#'    \item **burnin**: Numeric for the number of burnin before the first group.
#'      This is a key burnin as later processes often avoid burnin to ensure
#'      dependence / clean changes.
#'  }
#' @param parameters List of lists indicating parameters for each data group.
#'  Each process can use different parameters, given below.
#'  \itemize{
#'    \item **All**: All processes must contain:
#'      \itemize{
#'        \item **N**: Numeric indicating the number of observations.
#'        \item **process**: The process which should be used to generate the
#'          data. Options are the name for each: \code{bmotion}, \code{bbridge},
#'          \code{kl}, \code{ou}, and \code{far1}.
#'      }
#'    \item **bmotion**: Brownian motion contains:
#'      \itemize{
#'        \item **sd**: Numeric for the standard deviation of the observations.
#'      }
#'    \item **bbridge**: Brownian bridge contains:
#'      \itemize{
#'        \item **sd**: Numeric for the standard deviation of the observations.
#'      }
#'    \item **bbridge2**: The second order Brownian bridge contains:
#'      \itemize{
#'        \item **sd**: Numeric for the standard deviation of the Brownian motion
#'          used to build the data.
#'      }
#'    \item **kl**: Karhunen-Loeve expansion contains:
#'      \itemize{
#'        \item **distribution**: Distribution of the errors. Options include:
#'          \code{binomial}, \code{exponential}, \code{laplace}, \code{normal},
#'          \code{t} (add dof argument), and \code{cauchy}.
#'        \item **eigenvalues**: Numerics for the eigenvalues of the given
#'          distribution (value for each in basis).
#'        \item **mean**: Numeric for the mean of the group.
#'        \item **dependence**: Strength of dependence between observation.
#'        \item **basis**: fda \code{basisfd} object.
#'        \item **dof**: (Optional) Numeric for the degrees of freedom if using
#'          a \code{t} distribution.
#'      }
#'    \item **ou**: Ornstein-Uhlenbeck process requires:
#'      \itemize{
#'        \item **dependence**: Strength of dependence between observation.
#'      }
#'    \item **far1**: Functional autoregressive process of order 1 contains:
#'      \itemize{
#'        \item **dependence**: Strength of dependence between observation.
#'        \item **sd**: Numeric for the standard deviation of the observations.
#'        \item **vary**: Boolean if the starting value each observation should
#'          be 0 or vary. It does this by dropping first value. If resolution is
#'          given as a number, it can adjust so that the length is the same. If
#'          resolution is a vector, the resolution will be one smaller.
#'      }
#'  }
#'
#' @return funts object for the generated data
#' @export
#'
#' @examples
#' result <- generate_data(
#'   general=list('resolution'=10,'burnin'=100),
#'   parameters =list('bmotion'=list('N'=100, 'process'='bmotion', 'sd'=1),
#'                    'bbridge'=list('N'=100, 'process'='bbridge', 'sd'=1),
#'                    'kl'=list('process'='kl', 'N'=100,
#'                              'distribution'='Normal',
#'                              'eigenvalues'=1/1:4,
#'                              'mean'=0, 'dependence'=0,
#'                              'basis'=fda::create.bspline.basis(),
#'                              'sd'=1),
#'                    'ou'=list('N'=100, 'process'='ou', 'dependence'=0 ) ,
#'                    'far1'=list('N'=100, 'process'='far1', 'dependence'=0,
#'                                'sd'=1,'vary'=FALSE) )
#'                  )
#' result <- generate_data(
#'   general=list('resolution'=5,'burnin'=100),
#'   parameters =list(list('N'=100, 'dependence'=0, 'process'='far1', 'sd'=1,'vary'=TRUE)) )
generate_data <- function(general, parameters){

  if(length(general$resolution)==1){
    resolution <- seq(0,1, length.out=general$resolution)
  } else{
    resolution <- general$resolution
  }
  poss_process <- c('bbridge','bbridge2','bmotion','kl','ou','far1')
  data <- matrix()

  # Go through each defined groups
  for(i in 1:length(parameters)){
    # Parameters
    process <- .verify_input(parameters[[i]]$process, poss_process)
    N_data <- parameters[[i]]$N
    if(i==1){
      N <- N_data + general$burnin
    } else{
      N <- N_data
    }

    # Generate Data
    data_tmp <- switch(process,
           bbridge = {
             results_tmp <- generate_brownian_bridge(N = N, v = resolution,
                                                     sd = parameters[[i]]$sd)
             results_tmp$data[, (N-N_data+1):N]
           },
           bmotion = {
             results_tmp <- generate_brownian_motion(N = N, v = resolution,
                                                     sd = parameters[[i]]$sd)
             results_tmp$data[, (N-N_data+1):N]
           },
           kl = {
             results <- generate_karhunen_loeve(N = N,
                                     eigenvalues = parameters[[i]]$eigenvalues,
                                     basis = parameters[[i]]$basis,
                                     means = parameters[[i]]$mean,
                                     distribution = parameters[[i]]$distribution,
                                     resolution = resolution,
                                     dependence = parameters[[i]]$dependence,
                                     burnin = N-N_data,
                                     silent = FALSE,
                                     prev_eps = parameters[[i]]$prev_eps,
                                     dof=parameters[[i]]$dof)
             if(i < length(parameters))
               parameters[[i+1]]$prev_eps <- results$prev_eps

             results$data[, (N-N_data+1):N]
           },
           ou = {
             results <- .generate_ou(resolution = resolution, N = N,
                                     rho = parameters[[i]]$dependence)

             results$data[, (N-N_data+1):N]
           },
           far1 = {
             if(length(general$resolution)==1 & parameters[[i]]$vary){
               resolution1 <- seq(0, 1, length.out=general$resolution+1)
             }else if(parameters[[i]]$vary){
               stop(paste0('A specified resolution is given;',
                           'however, with paramter `vary=TRUE`, we must lose',
                           'a resolution point. It is unclear how to do so.',
                           'Specify resolution as a number or do not vary'))
             } else{
               resolution1 <- resolution
             }

             results <- generate_far1(N = N, resolution = resolution1,
                                       sd=parameters[[i]]$sd,
                                       dependence = parameters[[i]]$dependence,
                                       dropFirst=parameters[[i]]$vary)

             results$data[, (N-N_data+1):N]
           }
    )

    if(i==1){
      data <- data_tmp
    }else{
      data <- cbind(data,data_tmp)
    }
  }

  funts(data,intraobs = resolution,labels = 1:ncol(data))
}
