## generate a symmetric Q matrix

symmetric.Q.matrix <- function(K) {

  # Create an empty matrix of size K x K
  Q <- matrix(0, nrow = K, ncol = K)

  # transition rate for each state
  transition.rate <- 1 / (K)

  for (i in 1:K) {
    for (j in 1:K) {
      if (i != j) {
        Q[i, j] = transition.rate
      }
    }
    Q[i, i] = -(transition.rate*(K-1))
  }

  # Return the generated transition matrix
  return(Q)
}



## this is taken directly from phangorn.
get_gamma_rates <- function(alpha, k) {
  if (k == 1) return(1)
  # Compute quantiles to divide the gamma distribution into k categories
  quants <- qgamma((1:(k - 1)) / k, shape = alpha, rate = alpha)
  diff(c(0, pgamma(quants * alpha, alpha + 1), 1)) * k
}


