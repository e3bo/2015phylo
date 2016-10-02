#' Calculate the negative log likelihood of a tree occurring according
#' to a birth-death process.
#'
#' @param l Matrix of rates at which a type i individual gives birth
#' to a type j individual
#' @param m Vector of death rates for each type
#' @param psi Vector of death rates for each type where sampling
#' occurs at time of death. The probability of sampling at death is psi / (psi + m)
#' @param freq Expected probability of each type at the root of the
#' tree
#' @param phylo Phylo object with tree for which likelihood will be
#' calculated
#' @param survival Boolean indicator of whether likelihood is divided
#' by likelihood of tree being observed
#' @param unknown_states Boolean indicator of whether the type of the
#' tree tips is unknown
#' @param rtol Passed through to deSolve differential equation solver
#' @param atol Passed through to deSolve differential equation solver
#' @param cutoff Time in past from most recent tip beyond which
#' sampling probability assumed zero
#'
#' @export
calc_bd_nll <- function (l, m, psi, freq, phylo, survival = FALSE,
                        unknown_states = FALSE, rtol = 1e-12, atol = 1e-12,
                        cutoff = 10 ^ 12){
    maxpar <- 100
    summary <- get_times(phylo)
    out <- 10 ^ 1000
    ntypes <- nrow(l)

    bad_arg <- any(c(min(l, m, psi, freq) < 0, max(l, m, psi) > maxpar,
                   max(freq) > 1, length(freq) != ntypes - 1, sum(freq) > 1,
                   length(m) != ntypes, length(psi) != ntypes,
                   ncol(l) != ntypes))
    if (is.na(bad_arg)){
        bad_arg <- TRUE
    }
    if (!bad_arg) {
        lik <- try(get_subtree_lik(phylo, 1, l, m, psi, summary, unknown_states,
                                   rtol, atol, cutoff))
        if (class(lik) != "try-error") {
            pinds <- seq(1, ntypes)
            p <- lik[pinds]
            ginds <- seq(ntypes + 1, 2 * ntypes)
            g <- lik[ginds]
            freq <- c(freq, 1 - sum(freq))
            out <- sum(g * freq)
            if (survival) {
                out <- out / (1 - sum(p * freq))
            }
            if (any(c(class(out) != "numeric", !is.finite(out)))) {
                out <- 10 ^ 1000
            }
        }
    } else {
        stop("Invalid parameters")
    }
    -log(out)
}

get_times <- function (tree) {
    nodes <- sort(unique(c(tree$edge)))
    ttype <- (1:length(nodes)) * 0
    times <- ttype
    ttype[tree$edge[1, 1]] <- 1
    for (j in (tree$edge[1, 1] + 1):length(nodes)) {
        ttype[j] <- 1
        temp <- which(tree$edge[, 2] == j)
        ancestor <- tree$edge[temp, 1]
        times[j] <- times[ancestor] + tree$edge.length[temp]
    }
    for (j in 1:(tree$edge[1, 1] - 1)) {
        temp <- which(tree$edge[, 2] == j)
        ancestor <- tree$edge[temp, 1]
        times[j] <- times[ancestor] + tree$edge.length[temp]
    }
    maxt <- max(times)
    times <- -times + maxt
    out <- cbind(times, ttype)
    out
}

get_subtree_lik <- function (phylo, rootedge, l, m, psi, summary,
                             unknown_states, rtol, atol, cutoff) {
    ntypes <- length(m)
    newroot <- phylo$edge[rootedge, 2]
    newtrees <- which(phylo$edge[, 1] == newroot)
    tyoung <- summary[phylo$edge[rootedge, 2]]
    told <- summary[phylo$edge[rootedge, 1]]
    if (length(newtrees) == 0) {
        if (unknown_states == FALSE && phylo$states[newroot] > 0) {
            state <- phylo$states[newroot]
            initpsi <- numeric(ntypes)
            initpsi[state] <- psi[state]
        }
        else {
            initpsi <- psi
        }
        init <- rep_len(1, ntypes)
        inity1 <- solve_lik_unsampled(init, l, m, psi, c(0, tyoung), rtol, atol)
        if (told < cutoff) {
            res <- solve_lik(init = c(inity1, initpsi), l, m, psi,
                             c(tyoung, told), rtol, atol)
        }
        else {
            inity2 <- solve_lik(init = c(inity1, initpsi), l, m, psi,
                                c(tyoung, cutoff), rtol, atol)
            m <- m + psi
            psi <- rep_len(0, ntypes)
            res <- solve_lik(init = inity2, l, m, psi, c(cutoff, told), rtol,
                             atol)
        }
    }
    else {
        likleft <- get_subtree_lik(phylo, newtrees[1], l, m, psi, summary,
                                   unknown_states, rtol, atol, cutoff)
        likright <- get_subtree_lik(phylo, newtrees[2], l, m, psi, summary,
                                    unknown_states, rtol, atol, cutoff)
        pinds <- seq_along(m)
        ginds <- pinds + length(pinds)
        gleft <- likleft[ginds]
        gright <- likright[ginds]
        pbirth <- likleft[pinds]
        outer_prd <- outer(gleft, gright)
        gbirth <- rowSums(outer_prd * l + t(outer_prd) * l) / 2
        init1 <- c(pbirth, gbirth)
        if (tyoung > cutoff) {
            m <- m + psi
            psi <- rep_len(0, ntypes)
        } else if (told > cutoff) {
            init1 <- solve_lik(init = init1, l, m, psi, c(tyoung, cutoff),
                               rtol, atol)
            tyoung <- cutoff
            m <- m + psi
            psi <- rep_len(0, ntypes)
        }
        res <- solve_lik(init = init1, l, m, psi, c(tyoung, told), rtol,
                         atol)
    }
    res
}

solve_lik <- function (init, l, m, psi, times, rtol, atol) {
    pinds <- seq_along(m)
    ginds <- pinds + length(pinds)
    ode <- function(times, y, p) {
        with(as.list(c(y, p)), {
            p <- y[pinds]
            g <- y[ginds]
            dp <- - (rowSums(l) + m + psi) * p + (l * p) %*% p +  m
            dg <- - (rowSums(l) + m + psi) * g + (l * p) %*% g + (l * g) %*% p
            list(c(dp, dg))
        })
    }
    out <- deSolve::lsoda(init, times, ode, c(l, m, psi), rtol = rtol,
                          atol = atol)[2, -1]
    out
}

solve_lik_unsampled <- function (init, l, m, psi, times, rtol, atol) {
    ode <- function(times, y, p) {
        with(as.list(c(y, p)), {
            dy <- - (rowSums(l) + m + psi) * y + (l * y) %*% y +  m
            list(dy)
        })
    }
    p <- list(l, m, psi)
    out <- deSolve::lsoda(init, times, ode, p, rtol = rtol, atol = atol)[2, -1]
    out
}

#' Simulate tree according to a multi-type birth-death process
#'
#' @export
sim_bd_proc <- function (n, l, m, psi, init = 1){
    if (!requireNamespace("TreeSim", quietly = TRUE)) {
        stop(paste("The TreeSim package is needed for this function to work.",
                   "Please install it."),
             call. = FALSE)
    }
    maxpar <- 100
    ntypes <- nrow(l)
    check <- any(c(min(l, m, psi) < 0, max(l, m, psi) > maxpar,
                   length(m) != ntypes, length(psi) != ntypes,
                   ncol(l) != ntypes, init < 1))
    if (!check && !is.na(check)){
        lambdavector <- l
        deathvector <- psi + m
        sampprobvector <- psi / (psi + m)
        TreeSim::sim.bdtypes.stt.taxa(n, lambdavector, deathvector,
                                      sampprobvector, init=init)
    } else {
        stop("Invalid arguments")
    }
}

#' Generate parameter map for linear model interface to birth-death nll
#'
#' @export
gen_param_map <- function(n){
    function(x, w){
        n <- n
        scale <- exp(w[1])
        effects <- w[-1]
        stopifnot(nrow(x) == n^2)
        stopifnot(ncol(x) == length(effects))
        eta <- exp(x %*% effects)
        eta <- eta / mean(eta) * scale
        rate_matrix <- matrix(eta, nrow=n, ncol=n)
        ret <- list()
        ret$l <- rate_matrix
        ret$m <- rep(1, n) / 2
        ret$psi <- rep(1, n) / 2
        ret$survival <- FALSE
        ret$frequency <- c(1, rep(0, n - 2))
        ret
    }
}

#' Calculate negative log likelihood of linear model parameterization
#' of birth-death process
#'
#' @export
calc_bd_lm_nll <- function(w, x, y, param_map){
    pars <- param_map(x, w)
    calc_bd_nll(l=pars$l, m=pars$m, psi=pars$psi, freq=pars$freq, phylo=y,
                survival=pars$survival)
}
