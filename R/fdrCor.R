#' fdrcorr
## this function is used to calculate the false discovery rate (FDR) based on the list of P values.

fdrCor <- function(test_list) {
  # Use the p.adjust function from the stats package in R to adjust the p-values
  # Use the test_list and apply the false discovery rate method.
  p.adjust(test_list, method = 'fdr')
}

