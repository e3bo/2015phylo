
# R port of functions in "Spectral methods in Matlab" by Lloyd Trefethen

cheb <- function(N) {
  ## Computes differention matrix based on grid of Chebysev points
  ##
  ## Args:
  ##   N: degree of Chebysev polynomial interpolated on grid
  ##
  ## Value:
  ##   list containing D, the differentiation matrix, and x, the grid
  ##   of points
  if (N == 0) {
    return(list(D = 0, x = 1))
  }
  x <- cos(pi * (0:N) / N)
  c <- matrix(c(2, rep(1, N-1), 2), N + 1, 1) * (-1)^(0:N)
  X <- matrix(x, nrow = length(x), ncol = N + 1)
  dX <- X - t(X)
  D <- c %*% t(1 / c) / (dX + diag(N + 1))
  D <- D - diag(rowSums(D))
  list(D = D, x = x)
}
