myintegrator2 <- function (init, l, m, psi, times, rtol, atol) {
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

repnum_types <- function (l11, l12, l21, l22, death1, death2) {
    L <- l11 - l22 - death1 + death2
    c <- sqrt(L * L + 4 * l12 * l21)
    f1 <- (c + L) / (c + L + 2 * l12)
    R0 <- f1 * (l11 + l12) / death1 + (1 - f1) * (l22 + l21) / death2
    R0
}

mylik <- function (par, phylo, fix = rbind(c(0, 0), c(0, 0)), sampfrac,
    survival = 0, supercrit = 0, unknown_states = FALSE, rtol = 1e-12,
    atol = 1e-12, freq = 0, cutoff = 10 ^ 12, ntypes=2){
    prpar <- FALSE
    maxpar <- 100
    partemp <- vector()
    k <- 1
    for (i in 1:8) {
        index <- which(i == fix[1, ])
        if (length(index) > 0) {
            if (fix[2, index] >= 0) {
                partemp <- c(partemp, fix[2, index])
            }
            else {
                temp <- -fix[2, index]
                if (temp == 0.4) {
                  partemp <- c(partemp, partemp[3] * partemp[2] / partemp[1])
                }
                else {
                  partemp <- c(partemp, partemp[temp] * fix[3, index])
                }
            }
        }
        else {
            partemp <- c(partemp, par[k])
            k <- k + 1
        }
    }
    death <- partemp[5:6]
    l <- partemp[1:4]
    psi <- death * sampfrac
    m <- death * (1 - sampfrac)
    summary <- get_times2(phylo)

    out <- 10 ^ 10
    temp <- 1
    R0temp <- try(repnum_types(l[1], l[2], l[3], l[4], death[1], death[2]))
    if (supercrit == 1 && class(R0temp) == "numeric" && R0temp < 1) {
        temp <- 0
    }
    if (supercrit == 1 && class(R0temp) == "try-error") {
        temp <- 0
    }
    check <- any(c(length(which(partemp == "NaN")) > 0,  min(l, psi) < 0,
                   m < 0, max(l, m, psi) > maxpar, temp == 0))
    if (check) {
        out <- 10 ^ 10
    }
    else {
        lik <- try(bdss_num_help(phylo, 1, l, m, psi, summary, unknown_states,
                                 rtol, atol, cutoff))
        if (class(lik) != "try-error") {
            lamb_mu <- l[1] - l[4] - (m[1] + psi[1]) + (m[2] + psi[2])
            c <- sqrt(lamb_mu ^ 2 + 4 * l[2] * l[3])
            f1 <- (c + lamb_mu) / (c + lamb_mu + 2 * l[2])
            if (freq > 0) {
                f1 <- freq
            }
            out <- (lik[3] * (f1) + lik[4] * (1 - f1))
            if (survival == 1) {
                out <- out / (1 - (f1 * lik[1] + (1 - f1) * lik[2]))
            }
            if (any(c(class(out) != "numeric", out == "NaN", out == "Inf"))) {
                out <- 10 ^ 10
            }
        }
    }
    if (out == 10 ^ 10) {
        out <- 10 ^ 1000
    }
    if (prpar == TRUE) {
        print(par)
    }
    -log(out)
}

get_times2 <- function (tree) {
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

bdss_num_help <- function (phylo, rootedge, l, m, psi, summary, unknown_states,
                           rtol, atol, cutoff) {
    newroot <- phylo$edge[rootedge, 2]
    newtrees <- which(phylo$edge[, 1] == newroot)
    tyoung <- summary[phylo$edge[rootedge, 2]]
    told <- summary[phylo$edge[rootedge, 1]]
    if (length(newtrees) == 0) {
        if (unknown_states == FALSE && phylo$states[newroot] > 0) {
            state <- phylo$states[newroot]
            initpsi <- c(0, 0)
            initpsi[state] <- psi[state]
        }
        else {
            initpsi <- c(psi[1], psi[2])
        }
        inity1 <- myintegrator2(c(1, 1), l, m, psi, c(0, tyoung), rtol, atol)
        if (told < cutoff) {
            res <- myintegrator(init = c(inity1, initpsi), l, m, psi,
                                c(tyoung, told), rtol, atol)
        }
        else {
            inity2 <- myintegrator(init = c(inity1, initpsi), l, m, psi,
                                   c(tyoung, cutoff), rtol, atol)
            m[2] <- m[2] + psi[2]
            psi[2] <- 0
            m[1] <- m[1] + psi[1]
            psi[1] <- 0
            res <- myintegrator(init = inity2, l, m, psi, c(cutoff, told), rtol,
                                atol)
        }
    }
    else {
        likleft <- bdss_num_help(phylo, newtrees[1], l, m, psi, summary,
                                 unknown_states, rtol, atol, cutoff)
        likright <- bdss_num_help(phylo, newtrees[2], l, m, psi, summary,
                                  unknown_states, rtol, atol, cutoff)
        res1 <- c(likleft[1], likleft[3] * likright[3] * l[1])
        res1[2] <- res1[2] + likleft[3] * likright[4] * l[2] / 2
        res1[2] <- res1[2] + likleft[4] * likright[3] * l[2] / 2
        res2 <- c(likleft[2], likleft[4] * likright[4] * l[4])
        res2[2] <- res2[2] + likleft[3] * likright[4] * l[3] / 2
        res2[2] <- res2[2] + likleft[4] * likright[3] * l[3] / 2
        init1 <- c(res1[1], res2[1], res1[2], res2[2])
        if (tyoung > cutoff) {
            m[2] <- m[2] + psi[2]
            psi[2] <- 0
            m[1] <- m[1] + psi[1]
            psi[1] <- 0
        }
        if (tyoung < cutoff && told > cutoff) {
            init1 <- myintegrator(init = init1, l, m, psi, c(tyoung, cutoff),
                                  rtol, atol)
            tyoung <- cutoff
            m[2] <- m[2] + psi[2]
            psi[2] <- 0
            m[1] <- m[1] + psi[1]
            psi[1] <- 0
        }
        res <- myintegrator(init = init1, l, m, psi, c(tyoung, told), rtol,
                            atol)
    }
    res
}

myintegrator <- function (init, l, m, psi, times, rtol, atol) {
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
