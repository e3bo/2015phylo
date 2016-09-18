#' Calculate the likelihood of a tree occurring according to a
#' birth-death process.
#'
#' @param l Matrix of rates at which a type i individual gives birth
#' to a type j individual
#' @param m Vector of death rates for each type
#' @param psi Vector of sampling rates for each type
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
calc_bdlik <- function (l, m, psi, freq, phylo, survival = FALSE,
                        unknown_states = FALSE, rtol = 1e-12, atol = 1e-12,
                        cutoff = 10 ^ 12){
    maxpar <- 100
    summary <- get_times(phylo)
    out <- 10 ^ 1000
    ntypes <- nrow(l)

    check <- any(c(min(l, m, psi, freq) < 0, max(l, m, psi) > maxpar,
                   max(freq) > 1, length(freq) != ntypes - 1, sum(freq) > 1,
                   length(m) != ntypes, length(psi) != ntypes,
                   ncol(l) != ntypes))
    l <- as.numeric(l)
    if (! check || is.na(check)) {
        lik <- try(get_subtree_lik(phylo, 1, l, m, psi, summary, unknown_states,
                                 rtol, atol, cutoff))
        if (class(lik) != "try-error") {
            lik1inds <- seq(1, ntypes)
            lik1 <- lik[lik1inds]
            lik2inds <- seq(ntypes + 1, 2 * ntypes)
            lik2 <- lik[lik2inds]
            freq <- c(freq, 1 - sum(freq))
            out <- sum(lik2 * freq)
            if (survival) {
                out <- out / (1 - sum(lik1 * freq))
            }
            if (any(c(class(out) != "numeric", !is.finite(out)))) {
                out <- 10 ^ 1000
            }
        }
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

get_subtree_lik <- function (phylo, rootedge, l, m, psi, summary, unknown_states,
                           rtol, atol, cutoff) {
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
        inity1 <- solve_lik2(init, l, m, psi, c(0, tyoung), rtol, atol)
        if (told < cutoff) {
            res <- solve_lik(init = c(inity1, initpsi), l, m, psi,
                                c(tyoung, told), rtol, atol)
        }
        else {
            inity2 <- solve_lik(init = c(inity1, initpsi), l, m, psi,
                                   c(tyoung, cutoff), rtol, atol)
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
        res1 <- c(likleft[1], likleft[3] * likright[3] * l[1])
        res1[2] <- res1[2] + likleft[3] * likright[4] * l[2] / 2
        res1[2] <- res1[2] + likleft[4] * likright[3] * l[2] / 2
        res2 <- c(likleft[2], likleft[4] * likright[4] * l[4])
        res2[2] <- res2[2] + likleft[3] * likright[4] * l[3] / 2
        res2[2] <- res2[2] + likleft[4] * likright[3] * l[3] / 2
        init1 <- c(res1[1], res2[1], res1[2], res2[2])
        if (tyoung > cutoff) {
            psi <- rep_len(0, ntypes)
        }
        if (tyoung < cutoff && told > cutoff) {
            init1 <- solve_lik(init = init1, l, m, psi, c(tyoung, cutoff),
                                  rtol, atol)
            tyoung <- cutoff
            psi <- rep_len(0, ntypes)
        }
        res <- solve_lik(init = init1, l, m, psi, c(tyoung, told), rtol,
                            atol)
    }
    res
}

solve_lik <- function (init, l, m, psi, times, rtol, atol) {
    ode <- function(times, y, p) {
        lambda11 <- p[1]
        lambda12 <- p[2]
        lambda21 <- p[3]
        lambda22 <- p[4]
        mu1 <- p[5]
        mu2 <- p[6]
        psi1 <- p[7]
        psi2 <- p[8]
        yd1 <- mu1 - (lambda11 + lambda12 + mu1 + psi1) * y[1]
        yd1 <- yd1 + lambda11 * y[1] * y[1] + lambda12 * y[1] * y[2]
        yd2 <- mu2 - (lambda21 + lambda22 + mu2 + psi2) * y[2]
        yd2 <- yd2 + lambda21 * y[1] * y[2] + lambda22 * y[2] * y[2]
        yd3 <- - (lambda11 + lambda12 + mu1 + psi1) * y[3]
        yd3 <- yd3 + 2 * lambda11 * y[1] * y[3] + lambda12 * y[1] * y[4]
        yd3 <- yd3 + lambda12 * y[2] * y[3]
        yd4 <- - (lambda22 + lambda21 + mu2 + psi2) * y[4]
        yd4 <- yd4 + 2 * lambda22 * y[2] * y[4] + lambda21 * y[2] * y[3]
        yd4 <- yd4 + lambda21 * y[1] * y[4]
        list(c(yd1, yd2, yd3, yd4))
    }
    out <- lsoda(init, times, ode, c(l, m, psi), rtol = rtol,
                 atol = atol)[2, 2:5]
    out
}

solve_lik2 <- function (init, l, m, psi, times, rtol, atol) {
    ode <- function(times, y, p) {
        lambda11 <- p[1]
        lambda12 <- p[2]
        lambda21 <- p[3]
        lambda22 <- p[4]
        mu1 <- p[5]
        mu2 <- p[6]
        psi1 <- p[7]
        psi2 <- p[8]
        yd1 <- sum(mu1 - (lambda11 + lambda12 + mu1 + psi1) * y[1],
                   lambda11 * y[1] * y[1] + lambda12 * y[1] * y[2])
        yd2 <- sum(mu2 - (lambda21 + lambda22 + mu2 + psi2) * y[2],
                   lambda21 * y[1] * y[2] + lambda22 * y[2] * y[2])
        list(c(yd1, yd2))
    }
    p <- c(l, m, psi)
    out <- lsoda(init, times, ode, p, rtol = rtol, atol = atol)[2, 2:3]
    out
}
