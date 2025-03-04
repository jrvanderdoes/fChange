
recursive_segmentation <- function(X,
                                   alpha = 0.05,
                                       silent = FALSE,
                                       addAmt = 0) {
  if(ncol(X)<=1) return(NA)

  # Look for a single change
  potential_cp <- mean_change(X)

  # No Change Point Detected
  if (potential_cp$pvalue>alpha) {
    return()
  }

  # Display progress
  if (!silent) {
    cat(paste0(
      "ChangePoint Detected (", 1 + addAmt, "-", addAmt + ncol(X), " at ",
      addAmt + potential_cp$change, "): Segment Data and Re-Search\n"
    ))
  }

  # Search Recursively
  return(c(
    recursive_segmentation(
      X = as.data.frame(X[, 1:potential_cp$change]),
      alpha = alpha,
      silent = silent,
      addAmt = addAmt),
    potential_cp$change + addAmt,
    recursive_segmentation(
      X = as.data.frame(X[, (potential_cp$change + 1):ncol(X)]),
      alpha = alpha,
      silent = silent,
      addAmt = addAmt + potential_cp$change)
  ))
}
