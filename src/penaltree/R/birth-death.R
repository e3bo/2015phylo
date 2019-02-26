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
                        maxsteps = 1e4, 
                        cutoff = 10 ^ 12){
    maxpar <- 100

    out <- 10 ^ 1000
    ntypes <- nrow(l)

    bad_arg <- any(c(min(l, m, psi, freq) < 0, max(l, m, psi) > maxpar,
                     max(freq) > 1, length(freq) != ntypes,
                     abs(sum(freq) - 1) > .Machine$double.eps,
                     length(m) != ntypes, length(psi) != ntypes,
                     ncol(l) != ntypes))
    if (is.na(bad_arg)){
        bad_arg <- TRUE
    }
    if (!bad_arg) {
      phylor <- addroot(phylo, 0)
      summary <- get_times(phylor)
        #rootid <- length(phylor$tip.label) + 1
        #rootedge <- which (phylor$edge[, 1] == rootid)
        lik <- try(get_subtree_lik(phylor, 1L, l, m, psi, summary, unknown_states,
                                   rtol, atol, cutoff, maxsteps))
        if (class(lik) != "try-error") {
            pinds <- seq(1, ntypes)
            p <- lik[pinds]
            ginds <- seq(ntypes + 1, 2 * ntypes)
            g <- lik[ginds]
            #freq <- c(freq, 1 - sum(freq))
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
    ret <- -log(out)
    if(!is.finite(ret)){
                                        #browser()
        740
    } else {
        ret
    }
}

get_times <- function (tree) {
    nodes <- sort(unique(c(tree$edge)))
    ttype <- (1:length(nodes)) * 0
    times <- ttype
    rootid <- length(tree$tip.label) + 1
    ttype[rootid] <- 1
    for (j in (rootid + seq_len(length(nodes) - rootid))) {
        ttype[j] <- 1
        temp <- which(tree$edge[, 2] == j)
        ancestor <- tree$edge[temp, 1]
        times[j] <- times[ancestor] + tree$edge.length[temp]
    }
    for (j in 1:(rootid - 1)) {
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
                             unknown_states, rtol, atol, cutoff, maxsteps) {
    ntypes <- length(m)
    newroot <- phylo$edge[rootedge, 2]
    newtrees <- which(phylo$edge[, 1] == newroot)
    tyoung <- summary[phylo$edge[rootedge, 2]]
    told <- summary[phylo$edge[rootedge, 1]]
    if (length(newtrees) == 0) {
        if (unknown_states == FALSE && phylo$states[newroot] > 0) {
            state <- phylo$states[newroot]
            initpsi <- numeric(ntypes) + 1e-8
            initpsi[state] <- psi[state]
        }
        else {
            initpsi <- psi
        }
        init <- rep_len(1, ntypes)
        inity1 <- solve_lik_unsampled(init, l, m, psi, c(0, tyoung), rtol, atol, maxsteps)
        if (told < cutoff) {
            res <- solve_lik(init = c(inity1, initpsi), l, m, psi,
                             c(tyoung, told), rtol, atol, maxsteps)
        }
        else {
            inity2 <- solve_lik(init = c(inity1, initpsi), l, m, psi,
                                c(tyoung, cutoff), rtol, atol, maxsteps)
            m <- m + psi
            psi <- rep_len(0, ntypes)
            res <- solve_lik(init = inity2, l, m, psi, c(cutoff, told), rtol,
                             atol, maxsteps)
        }
    }
    else {
        likleft <- get_subtree_lik(phylo, newtrees[1], l, m, psi, summary,
                                   unknown_states, rtol, atol, cutoff, maxsteps)
        likright <- get_subtree_lik(phylo, newtrees[2], l, m, psi, summary,
                                    unknown_states, rtol, atol, cutoff, maxsteps)
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
                               rtol, atol, maxsteps)
            tyoung <- cutoff
            m <- m + psi
            psi <- rep_len(0, ntypes)
        }
        res <- solve_lik(init = init1, l, m, psi, c(tyoung, told), rtol,
                         atol, maxsteps)
    }
    res
}

solve_lik <- function (init, l, m, psi, times, rtol, atol, maxsteps) {
    pinds <- seq_along(m)
    ginds <- pinds + length(pinds)
    ode <- function(times, y, p) {
        with(as.list(c(y, p)), {
            pexp <- exp(y[pinds])
            gexp <- exp(y[ginds])
            dexpp <- - (rowSums(l) + m + psi) * pexp + (l * pexp) %*% pexp +  m
            dexpg <- - (rowSums(l) + m + psi) * gexp + (l * pexp) %*% gexp + (l * gexp) %*% pexp
            dp <- dexpp / pexp
            dg <- dexpg / gexp
            list(c(dp, dg))
        })
      }
    if (isTRUE(all.equal(times[2], times[1]))){
      return(init)
    } else {
      
      out <- try(deSolve::lsoda(log(init), times, ode, c(l, m, psi), rtol = rtol,
                                atol = atol, maxsteps = maxsteps)[2, -1])
    }
    if (inherits(out, "try-error")){
      browser()
    } else {
      exp(out)
    }
    return(exp(out))
}

solve_lik_unsampled <- function (init, l, m, psi, times, rtol, atol, maxsteps) {
    ode <- function(times, y, p) {
        with(as.list(c(y, p)), {
             yexp <- exp(y)
            dyexp <- - (rowSums(l) + m + psi) * yexp + (l * yexp) %*% yexp +  m
            list(dyexp / yexp)
        })
    }
    p <- list(l, m, psi)
    out <- try(deSolve::lsoda(log(init), times, ode, p, rtol = rtol, atol = atol, maxsteps = maxsteps)[2, -1])
    if (inherits(out, "try-error")){
         browser()
    } else {
         exp(out)
    }
}

#' Simulate tree according to a multi-type birth-death process
#'
#' @param n number of tips in simulated tree
#' @param l matrix of birth rates of each type for each type in process
#' @param m vector of death rates for each type. These deaths are not sampled.
#' @param psi vector of death rates for each type. These deaths are
#' followed by sampling.
#' @param init type of individual that starts the population
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
#' @param n number of types in birth-death model
#' @param ntrees number of trees that will be fitted to birth-death model
#' @param psampled proportion of deaths that are assumed to be
#' observed in birth-death model
#'
#' @export
gen_param_map <- function(n, ntrees, psampled=0.01){
    function(x, w){
        ret <- list(frequency=list())
        ret$m <- rep(exp(w[1]), n) * (1 - psampled)
        ret$psi <- rep(exp(w[1]), n) * psampled
        scale <- exp(w[2])
        start <- 3
        for (i in seq_len(ntrees)){
            inds <- seq(start, start - 1  + n - 1)
            foo <- c(1, exp(w[inds]))
            ret$frequency[[i]] <- foo / sum(foo)
            start <- start + n - 1
        }
        effects <- w[-seq(1, start - 1)]
        stopifnot(nrow(x) == n^2)
        stopifnot(ncol(x) == length(effects))
        eta <- exp(x %*% effects)
        eta <- eta / mean(eta) * scale
        rate_matrix <- matrix(eta, nrow=n, ncol=n)
        ret$l <- rate_matrix
        ret$survival <- FALSE
        ret
    }
}

#' Calculate negative log likelihood of linear model parameterization
#' of birth-death process
#'
#' @param w vector of parameters
#' @param x matrix of predictors
#' @param y tree assumed to come from birth-death process
#' @param param_map function to map w to list with parameters for calc_bd_nll()
#'
#' @export
calc_bd_lm_nll <- function(w, x, y, param_map){
    pars <- param_map(x, w)
    tmpf <- function(p1, p2) {
        calc_bd_nll(l=pars$l, m=pars$m, psi=pars$psi, freq=p1, phylo=p2,
                    survival=pars$survival)
    }
    nll <- mapply(tmpf, p1=pars$frequency, p2=y)
    sum(nll)
}
