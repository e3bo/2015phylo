
solve_lik_unsampled <- function (init, l, m, psi, times, rtol, atol) {
    ode <- function(times, y, p) {
        with(as.list(c(y, p)), {
            dy <- - (rowSums(l) + m + psi) * y + (l * y) %*% y +  m
            list(dy)
        })
    }
    p <- list(l, m, psi)
    try(out <- deSolve::lsoda(init, times, ode, p, rtol = rtol, atol = atol)[2, -1])
    if (inherits(out, "try-error")){
         browser()
    } else {
         out
    }
}

init <- c(1)
l <- matrix(0, ncol=1)
m <- 0.1
psi <- 0.1
rtol <- atol <- 1e-12
times <- c(0, 2)
ans1 <- solve_lik_unsampled(init = init, l = l, m = m, psi = psi, times = times, rtol = rtol, atol = atol)

layered_lik <- function(init, l, m, psi, times) {
    times <-  sort(times)[-1] - min(times)
    a <- sum(l) + m + psi
    C0 <- m / a    
    C1 <- 1 - C0
    p <- C1 * exp(-a * times) + C0
    names(p) <- times

    F <- function(t) C1 ^ 2 * exp(-a * t * 3) + 2 * C0 * C1 * exp(-a * t * 2) + C0 ^ 2 * exp(-a * t)
    p2 <- integrate(F, 1, times)$value * l[1,1] * exp(a * times)
    
    p2b <-  C1 ^ 2 * exp(-a * times * 3) / (-3 * a) + 2 * C0 * C1 * exp(-a * times * 2) / (-a * 2) + C0 ^ 2 * exp(-a * times) / (-a) + C1 ^ 2 * exp(-a * 3)/(3 * a) + 2 * C0 * C1 * exp(-a * 2) / (2 * a) + C0 ^ 2 * exp(-a ) / a
    p2b <- p2b * l[1,1] * exp(a * times)
    print(p2b)
    stopifnot(isTRUE(all.equal(p2, p2b)))
    p + p2b
}

ans2 <- layered_lik(init = init, l = l, m = m, psi = psi, times = times)
all.equal(ans1, ans2, check.attributes = FALSE)

l <- matrix(.001, ncol=1)
ans1 <- solve_lik_unsampled(init = init, l = l, m = m, psi = psi, times = times, rtol = rtol, atol = atol)
ans2 <- layered_lik(init = init, l = l, m = m, psi = psi, times = times)
all.equal(ans1, ans2, check.attributes = FALSE)
