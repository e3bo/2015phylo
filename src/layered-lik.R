
solve_lik_unsampled <- function (init, l, m, psi, times, rtol, atol) {
    ode <- function(times, y, p) {
        with(as.list(c(y, p)), {
            dy <- - (rowSums(l) + m + psi) * y + (l * y) %*% y +  m
            list(dy)
        })
    }
    p <- list(l, m, psi)
    deSolve::lsoda(init, times, ode, p, rtol = rtol, atol = atol)
}

init <- c(1)
l <- matrix(0, ncol=1)
m <- 0.2
psi <- 0.1
rtol <- atol <- 1e-12
times <- c(0, 3)
(ans1 <- solve_lik_unsampled(init = init, l = l, m = m, psi = psi, times = times, rtol = rtol, atol = atol))

d <- nrow(l)
dt <- 1e-4
tau <- max(times)
t0 <- 0
mesh <- seq(0, tau, by = dt)
npoints <- length(mesh)

B <- 10

phis <- list()

phit0_func1  <- function(t) diag(exp(-t * (rowSums(l) + m + psi)), nrow = d)
phit0 <- lapply(mesh, phit0_func1)
phit0inverse <- lapply(phit0, solve)
phitau <- lapply(phit0inverse, function(x) phit0[[npoints]] %*% x)
input_integrand <- sapply(phitau, function(x) x %*%  m)
p0 <- phit0[[npoints]] %*% init + sum(input_integrand * dt)
all.equal(ans1, as.numeric(p0), check.attributes = FALSE)

r <- (rowSums(l) + m + psi)
tau <- max(mesh)
expM0 <- diag(exp(-tau * r), rnow = d)
p0exact <- diag(exp(-tau * r), nrow = d) %*% init + diag((1 - exp(-r * tau)) / r, nrow= d) %*% m
all.equal(ans1, as.numeric(p0exact), check.attributes = FALSE)

p0exactint <- diag((1 - exp(-r * tau)) / r, nrow = d) %*% init + diag(tau / r - (1 - exp(-r * tau)) / r^2, nrow = d) %*% m





x <- solve_lik_unsampled(init = init, l = l, m = m, psi = psi, times = mesh, rtol = rtol, atol = atol)

all.equal(sum(x[,2] * dt), as.numeric(p0exactint))

expM1 <- diag(exp(l[1,1] * p0exactint))



phit0inverseex

input <- phitau

phitau <- phit0[length(mesh)]



for (i in seq(2, B + 1)){




}



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
