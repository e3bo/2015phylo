set.seed(41)

context("birth-death process")

test_that("Birth-death likelihood is consistent with TreePar function", {
    skip_if_not_installed("TreePar")
    skip_if_not_installed("TreeSim")

    lambda11 <- 15
    lambda12 <- 3
    lambda21 <- 1
    lambda22 <- 3
    death1 <- 4
    death2 <- 4
    sampprob1 <- 0.05
    sampprob2 <- 0.05
    l <- rbind(c(lambda11, lambda12), c(lambda21, lambda22))
    d <- c(death1, death2)
    s <- c(sampprob1, sampprob2)
    n <- 20
    init <- -1
    tree <- TreeSim::sim.bdtypes.stt.taxa(n, l, d, s, init)
    tree <- TreePar::addroot(tree, tree$root.edge)
    par <- c(2, 2, 2, 2)
    fix <- rbind(c(1, 6, 7, 8), c(15, -5, 0, 0), c(1, 1, 1, 1))
    tplik <- TreePar::LikTypesSTT(par=par, phylo=tree, fix=fix, sampfrac=s,
                                  survival=0, freq=0.1, posR=0)
    l <- rbind(c(15, 2), c(2, 2))
    m <- c(2, 2) * 0.95
    psi <- c(2, 2) * 0.05
    ptlik <- calc_bd_nll(l=l, m=m, psi=psi, freq=0.1, phylo=tree, survival=FALSE)
    expect_equal(unname(tplik), ptlik)
})

test_that("Birth-death likelihood runs with >2 types", {
    skip_if_not_installed("TreeSim")

    lambda11 <- 15
    lambda12 <- 3
    lambda21 <- 1
    lambda22 <- 3
    death1 <- 4
    death2 <- 4
    sampprob1 <- 0.05
    sampprob2 <- 0.05
    l <- rbind(c(lambda11, lambda12), c(lambda21, lambda22))
    d <- c(death1, death2)
    s <- c(sampprob1, sampprob2)
    n <- 20
    init <- -1
    tree <- TreeSim::sim.bdtypes.stt.taxa(n, l, d, s, init)
    tree <- TreePar::addroot(tree, tree$root.edge)

    l <- rbind(c(15, 2, 1), c(2, 2, 1), rep(1, 3))
    m <- c(2, 2, 2) * 0.95
    psi <- c(2, 2, 2) * 0.05
    ptlik <- calc_bd_nll(l=l, m=m, psi=psi, freq=c(0.9, 0.1), phylo=tree,
                        survival=FALSE)
    expect_equal(-36.2374678051509, ptlik)

    ntypes <- 14
    l <- matrix(2, ncol=ntypes, nrow=ntypes)
    m <- rep(1, ntypes)
    psi <- rep(0.05, ntypes)
    ptlik <- calc_bd_nll(l=l, m=m, psi=psi, freq=rep(1/ntypes, ntypes - 1),
                        phylo=tree, survival=FALSE)
    expect_equal(40.4926005495144, ptlik)
})

test_that("Score function has mean zero at the true parameter value", {
    skip_if_not_installed("TreeSim")
    skip_if_not_installed("TreePar")
    skip_if_not_installed("numDeriv")
    skip_if_not_installed("fizzlipuzzli")
    skip_on_cran()
    set.seed(1)
    l <- rbind(c(15, 3), c(1, 3))
    m <- c(1, 1) / 2
    psi <- c(1, 1) / 2
    capture.output(trees <- replicate(40, sim_bd_proc(n=40, l=l, m=m, psi=psi,
                                          init=1), simplify=FALSE))
    tmpf <- function(x) TreePar::addroot(x, x$root.edge)
    trees <- lapply(trees, tmpf)
    likwrap <- function(x, phylo){
        l[1, 1] <- x[1]
        l[2, 1] <- x[2]
        l[1, 2] <- x[3]
        l[2, 2] <- x[4]
        m[1] <- x[5]
        m[2] <- x[6]
        psi[1] <- x[7]
        psi[2] <- x[8]
        calc_bd_nll(l=l, m=m, psi=psi, freq=c(1), phylo=phylo, survival=FALSE)
    }
    get_score <- function(phylo){
        numDeriv::grad(likwrap, x=c(as.numeric(l), m, psi), phylo=phylo)
    }
    scores <- sapply(trees, get_score)
    htests <- apply(scores, 1, stats::t.test)
    ps <- sapply(htests, "[[", "p.value")
    expect_true(all(stats::p.adjust(ps) > 0.05))
})

test_that("Score function has mean zero for regression model", {
    skip_if_not_installed("TreeSim")
    skip_if_not_installed("TreePar")
    skip_if_not_installed("numDeriv")
    skip_on_cran()
    set.seed(1)
    l <- rbind(c(15, 3), c(1, 3))
    m <- c(1, 1) / 2
    psi <- c(1, 1) / 2
    capture.output(trees <- replicate(1, sim_bd_proc(n=40, l=l, m=m, psi=psi,
                                      init=1), simplify=FALSE))

    addroot <- function(x) TreePar::addroot(x, x$root.edge)
    trees <- lapply(trees, addroot)

    w <- c(log(mean(l)), 1)
    x <- matrix(as.numeric(log(l / mean(l))), nrow=4)

    pm <- gen_param_map(2)
    lm_nll <- calc_bd_lm_nll(w=w, x=x, y=trees[[1]], xw2pars=pm$xw2pars)
    nll <- calc_bd_nll(l=l, m=m, psi=psi, freq=c(1), phylo=trees[[1]], survival=FALSE)
    expect_equal(nll, lm_nll)

    load("testdata.rda")
    x <- x[seq(1, 169), ]
    x1 <- x[, 1, drop=FALSE]

    pm <- gen_param_map(13)
    w1 <- c(log(2), 0.5)
    pars <- pm$xw2pars(x=x1, w=w1)

    capture.output(trees <- replicate(1, sim_bd_proc(n=40, l=pars$l, m=pars$m, psi=pars$psi,
                                      init=1), simplify=FALSE))
    trees <- lapply(trees, addroot)
    lm_nll <- calc_bd_lm_nll(w=w1, x=x1, y=trees[[1]], xw2pars=pm$xw2pars)
    nll <- calc_bd_nll(l=pars$l, m=pars$m, psi=pars$psi, freq=pars$freq, phylo=trees[[1]], survival=FALSE)
    expect_equal(nll, lm_nll)
})
