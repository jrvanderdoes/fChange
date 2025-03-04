Standardize use of data, X, x ...

trim generic function

ccf
linear model?
  decompose
spectral analysis
- periodogram
expoentnail smoothing??
  decompose



format.percent <- function(x, ...) {
  ret <- formatC(vctrs::vec_data(x) * 100, digits = 1, format = "f")
  ret[is.na(x)] <- NA
  ret[!is.na(x)] <- paste0(ret[!is.na(x)], "%")
  return(ret)
}

print.fts

dimnames

`[`

`[.factor2` <- function(x, i) {
  new_factor2(NextMethod(), levels = attr(x, "levels"))
}

print.class <- function(x, ...) cat(format(x, ...), "\n")
