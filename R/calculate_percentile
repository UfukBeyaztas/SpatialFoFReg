calculate_percentile <- function(mat_list, percentile) {
  n <- nrow(mat_list[[1]])
  p <- ncol(mat_list[[1]])
  
  # Initialize result matrix
  percentile_matrix <- matrix(0, n, p)
  
  for (i in 1:n) {
    # Extract the i-th row from each matrix and combine into a matrix (bootstrap samples)
    rows_combined <- do.call(rbind, lapply(mat_list, function(mat) mat[i, ]))
    
    # Compute column-wise percentiles across bootstrap samples
    percentile_matrix[i, ] <- apply(
      rows_combined, 2,
      function(col) quantile(col, probs = percentile / 100)
    )
  }
  
  return(percentile_matrix)
}
