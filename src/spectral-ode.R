
# R port of functions in "Spectral methods in Matlab" by Lloyd Trefethen

cheb <- function(N, limits = c(-1, 1)) {
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
  scale <- max(limits) - min(limits)
  x <- (cos(pi * (0:N) / N) + (min(scale) - 1) )* scale / 2
  c <- matrix(c(2, rep(1, N-1), 2), N + 1, 1) * (-1)^(0:N)
  X <- matrix(x, nrow = length(x), ncol = N + 1)
  dX <- X - t(X)
  D <- c %*% t(1 / c) / (dX + diag(N + 1))
  D <- D - diag(rowSums(D))
  list(D = D, x = x)
}

p13 <- function() {
  ## solve linear BVP u_xx = exp(4x), u(-1)=u(1)=0
  ##
  ## Also plots the polynomial interpolation and the maximum difference with the exact solution.
  ##
  ## Value:
  ##   vector containing u, the value of the soluation on a Chebysev grid of 17 points

  N <- 16
  chebN <- cheb(N)
  D2 <- chebN$D %*% chebN$D
  D2 <- D2[2:N, 2:N]
  f <- exp(4 * chebN$x[2:N])
  u <- solve(D2, f)
  u <- c(0, u, 0)

  plot(chebN$x, u, pch = 16, ylab = "", xlab = "", xaxs = "i", yaxs = "i",
       ylim = c(-2.5, 0.5), tcl = 0.5, xpd = TRUE, col = "darkblue")
  grid()
  axis(3, tcl = 0.5, labels = FALSE)
  axis(4, tcl = 0.5, labels = FALSE)
  z <- poly(chebN$x, N)
  pf <- lm(u~ z)
  xx <-  seq(-1, 1, 0.01)
  zz <- predict(z, xx)
  uu <- cbind(1, zz) %*% coef(pf)
  lines(xx, uu, col = "darkblue")
  exact <- (exp(4*xx) - sinh(4)*xx - cosh(4) )/16
  lin <- signif(norm(uu - exact, type = "I"), 4)
  mtext(paste0("max err = ", lin))
  u
}

myp1 <- function() {
  ## solve linear BVP u_x = -2 x + .2, u(0)=1
  ##
  ## Also plots the polynomial interpolation and the maximum difference with the exact solution.
  ##
  ## Value:
  ##   vector containing u, the value of the soluation on a Chebysev grid of 31 points

  N <- 30
  chebN <- cheb(N)
  D <- chebN$D[1:N, 1:N] + diag(2, nrow = N)
  chebN$x <- (chebN$x + 1) * 5
  f <- rep(0.2, N)
  u <- solve(D, f)
  u <- c(u, 1)

  plot(chebN$x, u, pch = 16, ylab = "", xlab = "", xaxs = "i", yaxs = "i",
       tcl = 0.5, xpd = TRUE, col = "darkblue")
  grid()
  axis(3, tcl = 0.5, labels = FALSE)
  axis(4, tcl = 0.5, labels = FALSE)
  z <- poly(chebN$x, N - 1)
  pf <- lm(u~ z)
  xx <-  seq(-1, 1, 0.01)
  zz <- predict(z, xx)
  uu <- cbind(1, zz) %*% coef(pf)
  lines(xx, uu, col = "darkblue")
#  exact <- (exp(4*xx) - sinh(4)*xx - cosh(4) )/16
#  lin <- signif(norm(uu - exact, type = "I"), 4)
#  mtext(paste0("max err = ", lin))
  u
}

