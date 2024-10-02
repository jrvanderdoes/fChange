
```{r model_electricity, eval=FALSE}
set.seed(123456)
electricity_fd <- fda::Data2fd(1:24, electricity)

# sum(fda::eval.fd(1:24,electricity_fd)-electricity)
#   eval.fd: eval.basis(1:24,electricity_fd$basis) %*% electricity_fd$coefs

# 3 gets 95%, 7 gets 99%
nPCs <- 7
electricity_fpca <- fda::pca.fd(electricity_fd, nharm = nPCs)
electricity_fpca_comp <- electricity_fpca$scores

## Forecast Each
ts_dat <- list()
comps <- list()
comps_resids <- list()
for (i in 1:nPCs) {
  ts_dat[[i]] <- stats::ts(electricity_fpca_comp[, i], freq = 7)
  # comps[[i]] <- forecast::ets(ts_dat[[i]])
  # comps_resids[[i]] <- resid(comps[[i]])
  #comps[[i]] <- forecast::auto.arima(ts_dat[[i]])
  comps[[i]] <- forecast::ets(ts_dat[[i]])
  comps_resids[[i]] <- stats::resid(comps[[i]])
}

electricity_fpca_forecast <- do.call(cbind, comps_resids)

# Revert Back to fd
#   Want: 24 x 365
#     forecast: 365 x nPCs
#    coefs: 26  x nPCs
#      coefs %*% comps: 26 x 365
#    Eval: 24 x 26
#       eval %*% orig: 24 x 365
orig_coefs <- electricity_fpca$harmonics$coefs %*% t(electricity_fpca_forecast)
eval_fd_vals <- fda::eval.basis(1:24, electricity_fd$basis) %*% orig_coefs
```


```{r binary_segmantation, eval=FALSE}

tmp <- paperPlot(eval_fd_vals,res[[1]])
plotly::save_image(tmp,
                   "C:/Users/jerem/Downloads/elec_resids_changes_proj.png",
                   width=1600, height=800)

tmp <- paperPlot(electricity,res[[1]])
plotly::save_image(tmp,
                   "C:/Users/jerem/Downloads/elec_changes_proj.png",
                   width=1600, height=800)

# tmp <- paperPlot(electricity,res1)
# plotly::save_image(tmp,
#                    "C:/Users/jerem/Downloads/elec_changes.png",
#                    width=1600, height=800)
```


epca <- pcaExploration(electricity_fd,order = 7)
fpca <- pca_fit(funts(electricity_fd),TVE = 0.99,model = 'ets')
set.seed(654321)
res1 <- binary_segmentation(X = epca$residuals,
                            statistic='Tn',
                            method='Sim',
                            h_function = function(X){ncol(X)^(1/3)})
set.seed(654321)
res2 <- binary_segmentation(X = fpca$errors,
                            statistic='Tn',
                            method='Sim',
                            h_function = function(X){ncol(X)^(1/3)})

plot_fd(epca$residuals, res1[[1]])



The above detect a change point, and occasionally the placement; however, if we have multiple change points then another method is necessary. We use binary segmentation in this case.
```{r plot_function, eval=FALSE}
paperPlot <- function(data,CPs){
  curve_points = 1:nrow(data)
  plot_title = NULL
  val_axis_title = ""
  res_axis_title = ""
  FD_axis_title = ""
  FDReps = seq(as.Date("2014/01/01"), by = "day", length.out = 365)
  eye = list(x = -1.0, y = -3, z = 0.5)
  aspectratio=list(x=2,y=1.25,z=1.25)
  showticklabels=FALSE

  number <- length(data[1, ])
  valRange <- c(min(data), max(data))

  plotData <- data.frame(
    "resolution" = NA,
    "FDRep" = NA,
    "Color" = NA,
    "Value" = NA
  )[-1, ]

  # Color Group to first CP
  for (j in 1:min(CPs)) {
    plotData <- rbind(
      plotData,
      data.frame(
        "resolution" = curve_points,
        "FDRep" = FDReps[j],
        "Color" = 1,
        "Value" = data[, j]
      )
    )
  }
  # Color Group from last CP
  for (j in (max(CPs) + 1):number) {
    plotData <- rbind(
      plotData,
      data.frame(
        "resolution" = curve_points,
        "FDRep" = FDReps[j],
        "Color" = length(CPs) + 1,
        "Value" = data[, j]
      )
    )
  }

  # Color Additional Groups
  if (length(CPs) > 1) {
    for (i in 2:length(CPs)) {
      for (j in (CPs[i - 1] + 1):CPs[i]) {
        plotData <- rbind(
          plotData,
          data.frame(
            "resolution" = curve_points,
            "FDRep" = FDReps[j],
            "Color" = i,
            "Value" = data[, j]
          )
        )
      }
    }
  }

  scene <- list(
    camera = list(eye = eye),
    aspectmode = "manual",
    aspectratio = aspectratio
  )

  # Get Colors
  tmpColors <- RColorBrewer::brewer.pal(min(9, max(3, length(CPs) + 1)), "Set1")
  if (length(CPs) > 9) {
    tmpColors <- rep(tmpColors, ceiling(c(length(CPs) + 1) / 9))[1:(length(CPs) + 1)]
  }

  tidyr::`%>%`(
    tidyr::`%>%`(
      tidyr::`%>%`(
        plotly::plot_ly(plotData,
                        x = ~ as.factor(FDRep), y = ~resolution, z = ~Value,
                        type = "scatter3d", mode = "lines",
                        split = ~ as.factor(FDRep),
                        color = ~ as.factor(Color),
                        colors = tmpColors
        ),
        plotly::layout(
          scene = list(
            yaxis = list(
              title = res_axis_title,
              showticklabels = showticklabels
            ),
            xaxis = list(
              title = FD_axis_title,
              showticklabels = TRUE,
              dtick = "M2",
              tickformat="%b",
              tickfont =list(size=20)
            ),
            zaxis = list(
              title = val_axis_title,
              showticklabels = showticklabels
            )
          )
        )
      ),
      plotly::layout(title = plot_title, scene = scene)
    ),
    plotly::layout(showlegend = FALSE)
  )
}

```
