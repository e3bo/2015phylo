
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
l <- matrix(1, ncol=1)
m <- 0.2
psi <- 0.1
rtol <- atol <- 1e-12
times <- c(0, 3)
(ans1 <- solve_lik_unsampled(init = init, l = l, m = m, psi = psi, times = times, rtol = rtol, atol = atol))


ana <- function(lambda, mu, psi, t){
  (sqrt(-4*lambda*mu + (psi +
lambda + mu)^2)*cosh((t*sqrt(-4*lambda*mu + (psi + lambda + mu)^2))/2)
- (psi + lambda - mu)*sinh((t*sqrt(-4*lambda*mu + (psi + lambda +
mu)^2))/2))/(sqrt(-4*lambda*mu + (psi + lambda +
mu)^2)*cosh((t*sqrt(-4*lambda*mu + (psi + lambda + mu)^2))/2) + (psi -
lambda + mu)*sinh((t*sqrt(-4*lambda*mu + (psi + lambda + mu)^2))/2)) }

ans2 <- ana(l[1,1], m, psi, times[2])

stopifnot(isTRUE(all.equal(ans1, ans2, check.attributes = FALSE)))




library(rSymPy)


get_approx <- function(nlayers = 1){
    Var("xi, t, C1, C0, a, l")
    p <- sympy("C1 * exp(a * t) + C0")
    get_tdf <- function(f1, f2) {
        sympy("((", p, "+", f1, ")*(", p, "+", f2, ")).subs(t, xi)")
    }
    tdf <- get_tdf("0", "0")
    for (i in seq(1, nlayers)) {
        p2 <- sympy("integrate(exp(-a * xi) * ", tdf, ", (xi, 1, t)) * l * exp(a * t)")
        tdf <- get_tdf(p2, p2)
    }
    list(p, p2)
}



layered_lik <- function(init, l, m, psi, times, app) {
    times <-  sort(times)[-1] - min(times)
    a <- -(sum(l) + m + psi)
    C0 <- -m / a
    C1 <- 1 - C0
    num_eval <- function(exprstr) {
        eval(parse(text = exprstr))
    }
    t <- times
    num_eval(app[[1]]) + num_eval(app[[2]])
}

ans2 <- layered_lik(init = init, l = l, m = m, psi = psi, times = times)
all.equal(ans1, ans2, check.attributes = FALSE)

l <- matrix(.01, ncol=1)
system.time(ans1 <- solve_lik_unsampled(init = init, l = l, m = m, psi = psi, times = times, rtol = rtol, atol = atol))

system.time(app <- get_approx(nlayers = 2))
system.time(ans2 <- layered_lik(init = init, l = l, m = m, psi = psi, times = times, app = app))

all.equal(ans1, ans2, check.attributes = FALSE)
