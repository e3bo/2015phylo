
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

library(rSymPy)

layered_lik <- function(init, l, m, psi, times) {
    times <-  sort(times)[-1] - min(times)
    a <- -(sum(l) + m + psi)
    C0 <- -m / a
    C1 <- 1 - C0
    p <- C1 * exp(a * times) + C0
    Var("xi, t, C1, C0, a, l")
    f <- Sym("integrate( (C1 * exp(a * xi) + C0)**2, (xi, 1, t)) * l * exp(-a * t)")
    p2 <- sympy("(", f, ").subs(l, ",  as.numeric(l), ").subs(a,", a, ").subs(C0,", C0, ").subs(C1,", C1, ").subs(t,", times, ").evalf()")
    p + as.numeric(p2)
}

ans2 <- layered_lik(init = init, l = l, m = m, psi = psi, times = times)
all.equal(ans1, ans2, check.attributes = FALSE)

l <- matrix(.001, ncol=1)
ans1 <- solve_lik_unsampled(init = init, l = l, m = m, psi = psi, times = times, rtol = rtol, atol = atol)
ans2 <- layered_lik(init = init, l = l, m = m, psi = psi, times = times)
all.equal(ans1, ans2, check.attributes = FALSE)
