symmetric_Q_matrix <- function(K) {
  # Create an empty matrix of size K x K
  Q <- matrix(0, nrow = K, ncol = K)
  transition.prob <- 1 / (K)
  
  for (i in 1:K) {
    for (j in 1:K) {
      if (i != j) {
        Q[i, j] = transition.prob
      }
    }
    Q[i, i] = -(transition.prob*(K-1))
  }
  
  # Return the generated transition matrix
  return(Q)
}