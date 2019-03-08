
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
all.equal(ans1[2,2], as.numeric(p0), check.attributes = FALSE, tol = 1e-4)

r <- (rowSums(l) + m + psi)
tau <- max(mesh)
expM0 <- diag(exp(-tau * r), nrow = d)
p0exact <- expM0 %*% init + diag((1 - exp(-r * tau)) / r, nrow= d) %*% m
all.equal(ans1[2,2], as.numeric(p0exact), check.attributes = FALSE)

p0exactint <- diag((1 - exp(-r * tau)) / r, nrow = d) %*% init + diag(tau / r - (1 - exp(-r * tau)) / r^2, nrow = d) %*% m
x <- solve_lik_unsampled(init = init, l = l, m = m, psi = psi, times = mesh, rtol = rtol, atol = atol)
all.equal(sum(x[,2] * dt), as.numeric(p0exactint), tol = 1e-4)


# now for 2 layers

solve_lik2 <- function (init, l, m, psi, times, rtol, atol) {
    ode <- function(times, y, p) {
      with(as.list(c(y, p)), {
            y0 <- y[1]
            y1 <- y[2]
            dy0 <- - (rowSums(l) + m + psi) * y0 + m
            dy1 <- - (rowSums(l) + m + psi) * y1 + (l * y0) %*% y1 + m
            list(c(dy0, dy1))
        })
    }
    p <- list(l, m, psi)
    deSolve::lsoda(init, times, ode, p, rtol = rtol, atol = atol)
}

l2 <- matrix(1, ncol = 1, nrow = 1)
r2 <- (rowSums(l2) + m + psi)
dt <- 1e-2
mesh <- seq(0, 10, by = dt)
npoints <- length(mesh)
tau <- max(mesh)

ans2 <- solve_lik2(init = rep(1, ncol(l) * 2), l = l2,
                   m = m, psi = psi, times = c(0, tau), rtol = rtol, atol = atol)

expM0 <- diag(exp(-tau * r2), nrow = d)
p0exact <- expM0 %*% init + diag((1 - exp(-r2 * tau)) / r2, nrow= d) %*% m
all.equal(ans2[2, "1"], as.numeric(p0exact), check.attributes = FALSE)

p0exactint <- diag((1 - exp(-r2 * tau)) / r2, nrow = d) %*% init + diag(tau / r2 - (1 - exp(-r2 * tau)) / r2 ^ 2, nrow = d) %*% m

funA1 <- function(tau) {
  expM0 <- diag(exp(-tau * r2), nrow = d)
  p0exactint <- diag((1 - exp(-r2 * tau)) / r2, nrow = d) %*% init + diag(tau / r2 - (1 - exp(-r2 * tau)) / r2 ^ 2, nrow = d) %*% m
  expM1 <- diag(exp(l2 * p0exactint))
  expA1 <- expM1 %*% expM0
  expA1
}

expA1t <- lapply(mesh, funA1)
expA1tau <- lapply(expA1t, function(x) expA1t[[npoints]] %*%  (1/ x))
input_integrand <- sum (sapply(expA1tau, function(x) x %*% m)) * dt
p1approx <- expA1t[[npoints]] %*% init + input_integrand

all.equal(as.numeric(p1approx), as.numeric(ans2[2,"2"]), tol = 1e-2)


# now for 3 layers

solve_lik3 <- function (init, l, m, psi, times, rtol, atol) {
    ode <- function(times, y, p) {
      with(as.list(c(y, p)), {
            y0 <- y[1]
            y1 <- y[2]
            y2 <- y[3]
            dy0 <- - (rowSums(l) + m + psi) * y0 + m
            dy1 <- - (rowSums(l) + m + psi) * y1 + (l * y0) %*% y1 + m
            dy2 <- - (rowSums(l) + m + psi) * y2 + (l * y1) %*% y2 + m
            list(c(dy0, dy1, dy2))
        })
    }
    p <- list(l, m, psi)
    deSolve::lsoda(init, times, ode, p, rtol = rtol, atol = atol)
}

l3 <- matrix(2, ncol = 1, nrow = 1)
r3 <- (rowSums(l3) + m + psi)
dt <- 1e-3
mesh <- seq(0, 2, by = dt)
npoints <- length(mesh)
tau <- max(mesh)

ans3 <- solve_lik3(init = rep(1, ncol(l3) * 3), l = l3,
                   m = m, psi = psi, times = mesh, rtol = rtol, atol = atol)
ansInf <- solve_lik_unsampled(init = 1, l = l3, m = m , psi =psi, times =c(0, tau), rtol = rtol, atol = atol)
## 3 layers has a relative difference of less than 1e-3 with the infinite layer calculation

expM0 <- diag(exp(-tau * r3), nrow = d)
p0exact <- expM0 %*% init + diag((1 - exp(-r3 * tau)) / r3, nrow= d) %*% m
all.equal(ans3[2, "1"], as.numeric(p0exact), check.attributes = FALSE)

funA1 <- function(tau) {
  expM0 <- diag(exp(-tau * r3), nrow = d)
  p0exactint <- diag((1 - exp(-r3 * tau)) / r3, nrow = d) %*% init + diag(tau / r3 - (1 - exp(-r3 * tau)) / r3 ^ 2, nrow = d) %*% m
  expM1 <- diag(exp(l3 * p0exactint))
  expA1 <- expM1 %*% expM0
  expA1
}

expA1t <- lapply(mesh, funA1)
input_integrand <- list()
for (i in 1:npoints) {
  expA1taui <- lapply(expA1t[1:i], function(x) expA1t[[i]] %*%  (1 / x))
  print(i)
  input_integrand[[i]] <- sum (sapply(expA1taui, function(x) x %*% m)) * dt
}

p1approxv <- Map(function(x, y) x %*% init + y, expA1t, input_integrand)

p1_approx_rmse <- sqrt(mean((unlist(p1approxv) -  ans3[, "2"])^2))

funA2 <- function(pointno) {
  expM0 <- diag(exp(-dt * (pointno - 1) * r3), nrow = d)
  p1approxint <- sum(unlist(p1approxv)[1:pointno]) * dt
  expM1 <- diag(exp(l3 * p1approxint))
  expA1 <- expM1 %*% expM0
  expA1
}

expA2t <- lapply(seq_along(mesh), funA2)

input_integrand2 <- list()
for (i in 1:npoints) {
  expA2taui <- lapply(expA2t[1:i], function(x) expA2t[[i]] %*%  (1 / x))
  print(i)
  input_integrand2[[i]] <- sum (sapply(expA2taui, function(x) x %*% m)) * dt
}

p2approxv <- Map(function(x, y) x %*% init + y, expA2t, input_integrand2)
p2rmse <- sqrt(mean((ans3[,"3"] - unlist(p2approxv))^2))




## scraps


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
