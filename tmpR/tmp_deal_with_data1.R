library(ggplot2)

# CPsVals=complete_binary_segmentation(electricity, test_statistic_function=compute_Mn,
#                              fn=compute_Mn,
#                              cutoff_function = generalized_resampling,
#                              trim_function = trim_function,
#                              alpha = 0.05,
#                              iters=1000,
#                              M=1000)
# changepoint_verification(CPsVals=CPsVals, data=electricity,
#                          test_statistic_function=compute_Mn,
#                          fn=compute_Mn,
#                          cutoff_function = generalized_resampling,
#                          trim_function = trim_function,
#                          alpha = 0.05,
#                          iters=1000,
#                          M=1000)
#####################
results_elec <- elbow_method(data = electricity, test_statistic_function=compute_Mn,
                             fn=compute_Mn,
                             cutoff_function = generalized_resampling,
                             trim_function = trim_function,
                             alpha = 0.05,
                             iters=1000,
                             M=5000)
print(results_elec[[2]])
use_cp <- changepoint_verification(CPsVals = results_elec[[1]]$CP[2:4],
                                   data = electricity,
                                   test_statistic_function=compute_Mn,
                                   fn=compute_Mn,
                                   cutoff_function = generalized_resampling,
                                   trim_function = trim_function,
                                   M=5000, iters=1000, alpha=0.05)
plot_fd(electricity,CPs=results_elec[[1]]$CP[2:4])
plot_fd(electricity,CPs=use_cp)

##
results_elec_no <- elbow_method(data=electricity,
                                test_statistic_function=compute_Mn,
                                fn=compute_Mn,
                                cutoff_function = generalized_resampling,
                                trim_function = function(dat, ...){0},
                                alpha = 0.05,
                                iters=1000,
                                M=5000)
print(results_elec_no[[2]])
use_cp1 <- changepoint_verification(CPsVals = results_elec_no[[1]]$CP[2:11],
                                    data = electricity,
                                    test_statistic_function=compute_Mn,
                                    fn=compute_Mn,
                                    cutoff_function = generalized_resampling,
                                    trim_function = function(dat, ...){0},
                                    M=5000, iters=1000, alpha=0.05)
plot_fd(electricity,CPs=results_elec_no[[1]]$CP[2:11])
plot_fd(electricity,CPs=use_cp1)

##############
#     Take out Mean Change then Distribution
##############
# Overdetected so verifications removes (nearly) all
#results_elec_mean <- complete_binary_segmentation(data = electricity,
#                                                  changepoint_function=mean_change)
results_elec_mean1 <- elbow_method(data = electricity,
                                   test_statistic_function=compute_mean_stat,
                                   fn=compute_mean_stat,
                                   cutoff_function = compute_mean_cutoff,
                                   trim_function = function(dat, ...){0},
                                   alpha = 0.05)
results_elec_mean1[[2]]
use_cp2 <- changepoint_verification(CPsVals = results_elec_mean1[[1]]$CP[2:7],
                                    data = electricity,
                                    test_statistic_function=compute_mean_stat,
                                    fn=compute_mean_stat,
                                    cutoff_function = compute_mean_cutoff,
                                    trim_function = function(dat, ...){0},
                                    alpha = 0.05)
plot_fd(electricity,CPs=use_cp2)

segs <- list()
use_cp2_all <- c(0,use_cp2,ncol(electricity))
for(i in 1:(length(use_cp2_all)-1)){
  segs[[i]] <- complete_binary_segmentation(
                      data=electricity[, (use_cp2_all[i]+1):use_cp2_all[i+1]],
                      test_statistic_function = compute_Mn,
                      fn=compute_Mn,
                      cutoff_function = generalized_resampling,
                      trim_function = function(dat, ...){0},
                      alpha = 0.05, iters=1000, M=5000,
                      final_verify = T,
                      errorType='Tr') + use_cp2_all[[i]]
}
plot_fd(electricity,CPs= c(use_cp2,unlist(segs)[!is.na(unlist(segs))]))
saveRDS(c(use_cp2,unlist(segs)[!is.na(unlist(segs))]),
        file='C:/Users/Jeremy/Downloads/elbowplots+/New/elbow_mean_v_full_Mn_v.RDS')

## Method 2
use_cp_elec_all <- c(0,use_cp2,ncol(electricity))
elec_demeaned <- electricity
for(i in 1:(length(use_cp_elec_all)-1)){
  section_mean <- rowMeans(electricity[,(use_cp_elec_all[i]+1):use_cp_elec_all[i+1]])

  elec_demeaned[,(use_cp_elec_all[i]+1):use_cp_elec_all[i+1]] <-
    electricity[,(use_cp_elec_all[i]+1):use_cp_elec_all[i+1]] - section_mean
}

dist_seg <- complete_binary_segmentation(
  data=elec_demeaned,
  test_statistic_function = compute_Mn,
  fn=compute_Mn,
  cutoff_function = generalized_resampling,
  trim_function = function(dat, ...){0},
  alpha = 0.05, iters=1000, M=5000,
  final_verify = T,
  errorType='Tr')

if(!is.na(dist_seg))
  plot_fd(electricity,CPs= c(use_cp2,dist_seg))
if(is.na(dist_seg))
  plot_fd(electricity,CPs= use_cp2)

##############
dist_seg <- complete_binary_segmentation(
  data=electricity,
  test_statistic_function = compute_Tn,
  fn=compute_Tn,
  cutoff_function = generalized_resampling,
  trim_function = function(dat, ...){0},
  alpha = 0.05, iters=1000, M=5000,
  final_verify = T,
  errorType='Tr')
dist_seg

####################################
results_syd <- elbow_method(aust_tmp2, test_statistic_function=compute_Mn,
                            fn=compute_Mn,
                            cutoff_function = generalized_resampling,
                            trim_function = trim_function,
                            alpha = 0.05,
                            iters=1000,
                            M=5000)
print(results_syd[[2]])
use_cp4 <- changepoint_verification(CPsVals = results_syd[[1]]$CP[2:3],
                                    data = aust_tmp2,
                                    test_statistic_function=compute_Mn,
                                    fn=compute_Mn,
                                    cutoff_function = generalized_resampling,
                                    trim_function = trim_function,
                                    M=5000, iters=1000, alpha=0.05)
plot_fd(aust_tmp2,CPs=results_syd[[1]]$CP[2:3])
plot_fd(aust_tmp2,CPs=use_cp4)


results_syd_no <- elbow_method(aust_tmp2, test_statistic_function=compute_Mn,
                               fn=compute_Mn,
                               cutoff_function = generalized_resampling,
                               trim_function = function(dat, ...){0},
                               alpha = 0.05,
                               iters=1000,
                               M=5000)
results_syd_no[[2]]
use_cp5 <- changepoint_verification(CPsVals = results_syd_no[[1]]$CP[2:5],
                                    data = aust_tmp2,
                                    test_statistic_function=compute_Mn,
                                    fn=compute_Mn,
                                    cutoff_function = generalized_resampling,
                                    trim_function = function(dat, ...){0},
                                    M=5000, iters=1000, alpha=0.05)
plot_fd(aust_tmp2,CPs=results_syd_no[[1]]$CP[2:5])
plot_fd(aust_tmp2,CPs=use_cp5)


##############
#     Take out Mean Change then Distribution
##############
results_syd_mean1 <- elbow_method(data = aust_tmp2,
                                   test_statistic_function=compute_mean_stat,
                                   fn=compute_mean_stat,
                                   cutoff_function = compute_mean_cutoff,
                                   trim_function = function(dat, ...){0},
                                   alpha = 0.05)
results_syd_mean1[[2]]
use_cp_syd <- changepoint_verification(CPsVals = results_syd_mean1[[1]]$CP[2:3],
                                    data = aust_tmp2,
                                    test_statistic_function=compute_mean_stat,
                                    fn=compute_mean_stat,
                                    cutoff_function = compute_mean_cutoff,
                                    trim_function = function(dat, ...){0},
                                    alpha = 0.05)


## Method 1
segs_syd <- list()
use_cp_syd_all <- c(0,use_cp_syd,ncol(aust_tmp2))
for(i in 1:(length(use_cp_syd_all)-1)){
  segs_syd[[i]] <- complete_binary_segmentation(
    data=aust_tmp2[, (use_cp_syd_all[i]+1):use_cp_syd_all[i+1]],
    test_statistic_function = compute_Mn,
    fn=compute_Mn,
    cutoff_function = generalized_resampling,
    trim_function = function(dat, ...){0},
    alpha = 0.05, iters=1000, M=5000,
    final_verify = T,
    errorType='Tr') + use_cp_syd_all[[i]]
}
plot_fd(aust_tmp2,CPs= c(use_cp_syd,unlist(segs_syd)[!is.na(unlist(segs_syd))]))
#saveRDS(c(use_cp_syd,unlist(segs_syd)[!is.na(unlist(segs_syd))]),
#        file='C:/Users/Jeremy/Downloads/elbowplots+/New/syd_elbow_mean_v_full_Mn_v.RDS')

## Method 2
use_cp_syd_all <- c(0,use_cp_syd,ncol(aust_tmp2))
aust_tmp2_demeaned <- aust_tmp2
#aust_tmp2_demeaned1 <- aust_tmp2
for(i in 1:(length(use_cp_syd_all)-1)){
  section_mean <- rowMeans(aust_tmp2[,(use_cp_syd_all[i]+1):use_cp_syd_all[i+1]])

  #aust_tmp2_demeaned[,(use_cp_syd_all[i]+1):use_cp_syd_all[i+1]] <-
  #  aust_tmp2[,(use_cp_syd_all[i]+1):use_cp_syd_all[i+1]] - section_mean


  aust_tmp2_demeaned[,(use_cp_syd_all[i]+1):use_cp_syd_all[i+1]] <-
    aust_tmp2[,(use_cp_syd_all[i]+1):use_cp_syd_all[i+1]] - mean(section_mean)
}

dist_seg <- complete_binary_segmentation(
    data=aust_tmp2_demeaned,
    test_statistic_function = compute_Mn,
    fn=compute_Mn,
    cutoff_function = generalized_resampling,
    trim_function = function(dat, ...){0},
    alpha = 0.05, iters=1000, M=5000,
    final_verify = T,
    errorType='Tr')
dist_seg
if(!is.na(dist_seg))
  plot_fd(aust_tmp2,CPs= c(use_cp_syd,dist_seg))
if(is.na(dist_seg))
  plot_fd(aust_tmp2,CPs= use_cp_syd)
