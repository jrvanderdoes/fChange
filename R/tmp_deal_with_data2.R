
data_KL <- generate_data_fd(ns = c(100,50,100),
                             eigsList = list(c(3,2,1,0.5),
                                             c(3,2,1,0.5)*10,
                                             c(3,2,1,0.5)),
                             basesList = list(
                               fda::create.bspline.basis(nbasis=4, norder=4),
                               fda::create.bspline.basis(nbasis=4, norder=4),
                               fda::create.bspline.basis(nbasis=4, norder=4)),
                             meansList = c(0,0.5,0),
                             distsArray = c('Normal','Normal','Normal'),
                             evals = seq(0,1,0.05),
                             kappasArray = c(0, 0, 0))

results_KL <- elbow_method(data=data_KL,
                            test_statistic_function=compute_Mn, fn=compute_Mn,
                            cutoff_function = generalized_resampling,
                            trim_function = trim_function,
                            alpha = 0.05, errorType = 'CE',
                            M=10000)
print(results_KL[[2]])
results_KL[[1]]$CP[2:3]

CPsVals_KL <- complete_binary_segmentation(data = data_KL,
                    test_statistic_function=compute_Mn,
                    fn=compute_Mn, cutoff_function = generalized_resampling,
                    trim_function = trim_function,
                    alpha = 0.05, iters=1000,  M=1000)
CPsVals_KL

##################################
results_elec <- elbow_method(data=electricity,
                           test_statistic_function=compute_Mn, fn=compute_Mn,
                           cutoff_function = generalized_resampling,
                           trim_function = trim_function,
                           alpha = 0.05, errorType = 'CE',
                           M=1000)
print(results_elec[[2]])
results_elec[[1]]$CP[2:3]

CPsVals_elec <- complete_binary_segmentation(data = electricity,
                   test_statistic_function=compute_Mn,
                   fn=compute_Mn, cutoff_function = generalized_resampling,
                   trim_function = trim_function,
                   alpha = 0.05, iters=2500,  M=2500)
CPsVals_elec
plot_fd(electricity, CPs=CPsVals_elec)
tmp <- complete_binary_segmentation(data = electricity[,183:ncol(electricity)],
                                              test_statistic_function=compute_Mn,
                                              fn=compute_Mn, cutoff_function = generalized_resampling,
                                              trim_function = trim_function,
                                              alpha = 0.05, iters=2000,  M=2000)
plot_fd(electricity, CPs=c(CPsVals_elec,CPsVals_elec+tmp))
