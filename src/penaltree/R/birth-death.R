#' Calculate the likelihood of a tree occurring according to a
#' birth-death process.
#'
#' @param l Matrix of rates at which a type i individual gives birth
#' to a type j individual
#' @param m Vector of death rates for each type
#' @param psi Vector of sampling probability for each type. Sampling
#' occurs at time of death with this probability
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
    if (! check || is.na(check)) {
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

sim_bd_process <- function (n, lambdavector, deathvector, sampprobvector, init = -1, 
    EI = FALSE, eliminate = 0) 
{
    muvector <- deathvector * (1 - sampprobvector)
    psivector <- deathvector * sampprobvector
    extincttree = 1
    if ((init == -1) && (length(deathvector) == 2)) {
        init <- 2
        lamb <- lambdavector[1, 1] - lambdavector[2, 2] - deathvector[1] + 
            deathvector[2]
        c <- sqrt(lamb^2 + 4 * lambdavector[1, 2] * lambdavector[2, 
            1])
        f1 <- (c + lamb)/(c + lamb + 2 * lambdavector[1, 2])
        r <- runif(1, 0, 1)
        if (r < f1) {
            init <- 1
        }
    }
    if ((init == -1) && (length(deathvector) != 2)) {
        init <- sample(1:length(deathvector), 1)
    }
    if (EI == TRUE) {
        init <- 1
    }
    while (extincttree == 1) {
        edge <- c(-1, -2)
        leaves <- c(-2)
        types <- c(init)
        sampled <- vector()
        typessampled <- vector()
        timecreation <- c(0, 0)
        extinct <- vector()
        time <- 0
        maxspecies <- -2
        edge.length <- c(0)
        extincttree = 0
        stop = 0
        while (stop == 0) {
            if (length(leaves) == 0) {
                phy2 = 0
                extincttree = 1
                print("extinct tree")
                stop = 1
            }
            else {
                sumrates <- vector()
                for (i in 1:length(lambdavector[, 1])) {
                  sumrates <- c(sumrates, (length(which(types == 
                    i)) * (sum(lambdavector[i, ]) + muvector[i] + 
                    psivector[i])))
                }
                timestep <- rexp(1, sum(sumrates))
                time = time + timestep
                r <- runif(1, 0, sum(sumrates))
                chosentype <- min(which(cumsum(sumrates) > r))
                species <- sample(leaves[which(types == chosentype)], 
                  1)
                lambda <- sum(lambdavector[chosentype, ])
                gamma <- 0
                if (EI == TRUE) {
                  if (chosentype == 1) {
                    gamma <- lambda
                    lambda <- 0
                  }
                }
                mu <- muvector[chosentype]
                psi <- psivector[chosentype]
                del <- which(leaves == species)
                specevent <- runif(1, 0, 1)
                edgespecevent <- which(edge == species) - length(edge.length)
                if ((lambda/(lambda + gamma + mu + psi)) > specevent) {
                  edge.length[edgespecevent] <- time - timecreation[-species]
                  edge <- rbind(edge, c(species, maxspecies - 
                    1))
                  edge <- rbind(edge, c(species, maxspecies - 
                    2))
                  edge.length <- c(edge.length, 0, 0)
                  r <- runif(1, 0, lambda)
                  newtype <- min(which(cumsum(lambdavector[chosentype, 
                    ]) > r))
                  leaves <- c(leaves, maxspecies - 1, maxspecies - 
                    2)
                  types <- c(types, chosentype, newtype)
                  maxspecies <- maxspecies - 2
                  leaves <- leaves[-del]
                  types <- types[-del]
                  timecreation <- c(timecreation, time, time)
                }
                else if (((lambda + gamma)/(lambda + gamma + 
                  mu + psi)) > specevent) {
                  types[del] <- 2
                }
                else if (((lambda + gamma + psi)/(lambda + gamma + 
                  mu + psi)) > specevent) {
                  sampled <- c(sampled, leaves[del])
                  if (EI == TRUE && length(typessampled) < eliminate) {
                    typessampled <- c(typessampled, 1)
                  }
                  else {
                    typessampled <- c(typessampled, chosentype)
                  }
                  leaves <- leaves[-del]
                  types <- types[-del]
                  edge.length[edgespecevent] <- time - timecreation[-species]
                  if (length(sampled) == n) {
                    stop = 1
                  }
                }
                else {
                  extinct <- c(extinct, leaves[del])
                  leaves <- leaves[-del]
                  types <- types[-del]
                  edge.length[edgespecevent] <- time - timecreation[-species]
                }
            }
        }
    }
    while (length(leaves) > 0) {
        del <- 1
        extinct <- c(extinct, leaves[del])
        k = which(edge == leaves[del]) - length(edge.length)
        edge.length[k] <- time - timecreation[-leaves[del]]
        leaves <- leaves[-del]
    }
    for (j in 1:length(extinct)) {
        del <- which(edge == extinct[j]) - length(edge.length)
        surpress <- edge[del, 1]
        edge.length <- edge.length[-del]
        edge <- edge[-del, ]
        del2 <- which(edge[, 1] == surpress)
        modify <- which(edge[, 2] == surpress)
        edge[modify, 2] <- edge[del2, 2]
        edge.length[modify] <- edge.length[modify] + edge.length[del2]
        edge.length <- edge.length[-del2]
        edge <- edge[-del2, ]
    }
    leaf = 1
    interior = length(sampled) + 1
    edgetemp <- edge
    typessampledNew <- vector()
    temp <- unique(c(edge))
    temp <- sort(temp, decreasing = TRUE)
    for (j in temp) {
        if (sum(match(sampled, j, 0)) == 0 || j == -1) {
            posvalue <- interior
            interior <- interior + 1
        }
        else {
            typessampledNew <- c(typessampledNew, typessampled[which(sampled == 
                j)])
            posvalue <- leaf
            leaf <- leaf + 1
        }
        replacel <- which(edge == j)
        if (length(replacel) > 0) {
            for (k in 1:length(replacel)) {
                if ((replacel[k] - 1) < length(edge.length)) {
                  edge[replacel[k], 1] <- posvalue
                }
                else {
                  edge[(replacel[k] - length(edge.length)), 2] <- posvalue
                }
            }
        }
    }
    phy <- list(edge = edge)
    phy$tip.label <- paste("t", sample(length(sampled)), sep = "")
    phy$edge.length <- edge.length
    phy$Nnode <- length(sampled)
    phy$states <- typessampledNew
    class(phy) <- "phylo"
    br <- sort(getx(phy, sersampling = 1), decreasing = TRUE)
    phy$root.edge <- br[1] - br[2]
    phy <- collapse.singles(phy)
    phy2 <- phy
    phy2
}
