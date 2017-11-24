
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
    C0 <- m / (sum(l) + m + psi)
    C1 <- 1 - C0
    p <- C1 * exp(-(sum(l) + m + psi) * times) + C0
    names(p) <- times
    p
}

ans2 <- layered_lik(init = init, l = l, m = m, psi = psi, times = times)
all.equal(ans1, ans2, check.attributes = FALSE)
