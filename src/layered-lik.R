
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
    num_eval <- function(expr) {
        ret <- sympy("(", expr, ").subs(l, ",  as.numeric(l), ").subs(a,", a, ").subs(C0,", C0, ").subs(C1,", C1, ").subs(t,", times, ").evalf()")
        as.numeric(ret)
    }
    Var("xi, t, C1, C0, a, l")
    p <- sympy("C1 * exp(a * t) + C0")
    get_tdf <- function(f1, f2) {
        sympy("((", f1, ")*(", f2, ")).subs(t, xi)")
    }
    tdf <- get_tdf(p, p)
    p2 <- Sym("integrate(exp(-a * xi) * ", tdf, ", (xi, 1, t)) * l * exp(a * t)")
    tdf <- get_tdf(p2, p)
    p3 <- Sym("integrate(exp(-a * xi) * ", tdf, ", (xi, 1, t)) * l * exp(a * t)")
    tdf <- get_tdf(p2, p2)
    p4 <- Sym("integrate(exp(-a * xi) * ", tdf, ", (xi, 1, t)) * l * exp(a * t)")
    num_eval(psym) + num_eval(p2) + num_eval(p3) * 2 + num_eval(p4)
}

ans2 <- layered_lik(init = init, l = l, m = m, psi = psi, times = times)
all.equal(ans1, ans2, check.attributes = FALSE)

l <- matrix(.01, ncol=1)
ans1 <- solve_lik_unsampled(init = init, l = l, m = m, psi = psi, times = times, rtol = rtol, atol = atol)
system.time(ans2 <- layered_lik(init = init, l = l, m = m, psi = psi, times = times))
all.equal(ans1, ans2, check.attributes = FALSE)
