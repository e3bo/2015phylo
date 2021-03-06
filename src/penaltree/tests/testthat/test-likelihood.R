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
    treer <- TreePar::addroot(tree, tree$root.edge)
    par <- c(2, 2, 2, 2)
    fix <- rbind(c(1, 6, 7, 8), c(15, -5, 0, 0), c(1, 1, 1, 1))
    tplik <- TreePar::LikTypesSTT(par=par, phylo=treer, fix=fix, sampfrac=s,
                                  survival=0, freq=0.1, posR=0)
    l <- rbind(c(15, 2), c(2, 2))
    m <- c(2, 2) * 0.95
    psi <- c(2, 2) * 0.05
    ptlik <- calc_bd_nll(l=l, m=m, psi=psi, freq=c(0.1, 0.9),
                         phylo=tree, survival=FALSE)
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

    l2 <- rbind(c(15, 2, 1), c(2, 2, 1), rep(1, 3))
    m2 <- c(2, 2, 2) * 0.95
    psi2 <- c(2, 2, 2) * 0.05
    ptlik <- calc_bd_nll(l=l2, m=m2, psi=psi2, freq=c(0.9, 0.1, 0), phylo=tree,
                        survival=FALSE)
    expect_equal(-36.2374678051509, ptlik)

    ntypes <- 14
    l <- matrix(2, ncol=ntypes, nrow=ntypes)
    m <- rep(1, ntypes)
    psi <- rep(0.05, ntypes)
    ptlik <- calc_bd_nll(l=l, m=m, psi=psi, freq=rep(1/ntypes, ntypes),
                        phylo=tree, survival=FALSE)
    expect_equal(40.4926005495144, ptlik)
})

test_that("Score function has mean zero at the true parameter value", {
    skip_if_not_installed("TreeSim")
    skip_if_not_installed("TreePar")
    skip_on_cran()
    set.seed(1)
    l <- rbind(c(15, 3), c(1, 3))
    m <- c(1, 1) / 2
    psi <- c(1, 1) / 2
    capture.output(trees <- replicate(20, sim_bd_proc(n=40, l=l, m=m, psi=psi,
                                          init=1), simplify=FALSE))
    likwrap <- function(x, phylo){
        l[1, 1] <- x[1]
        l[2, 1] <- x[2]
        l[1, 2] <- x[3]
        l[2, 2] <- x[4]
        m[1] <- x[5]
        m[2] <- x[6]
        psi[1] <- x[7]
        psi[2] <- x[8]
        calc_bd_nll(l=l, m=m, psi=psi, freq=c(1, 0), phylo=phylo, survival=FALSE)
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
    skip_on_cran()
    set.seed(1)
    l <- rbind(c(15, 3), c(1, 3))
    m <- c(1, 1) / 2
    psi <- c(1, 1) / 2
    capture.output(trees <- replicate(1, sim_bd_proc(n=40, l=l, m=m, psi=psi,
                                      init=1), simplify=FALSE))

    w <- c(log(m[1] + psi[1]), log(mean(l)), 0, 1)
    x <- matrix(as.numeric(log(l / mean(l))), nrow=4)

    pm <- gen_param_map(n = 2, ntrees = 1, psampled = 0.5)
    lm_nll <- calc_bd_lm_nll(w=w, x=x, y=trees[1], param_map=pm)
    nll <- calc_bd_nll(l=l, m=m, psi=psi, freq=c(0.5, 0.5), phylo=trees[[1]],
                       survival=FALSE)
    expect_equal(nll, lm_nll)

    load("testdata.rda")
    x <- x[seq(1, 169), ]
    x1 <- x[, c(1, 2), drop=FALSE]

    pm <- gen_param_map(13, ntrees=1, psampled=0.5)
    w1 <- c(log(m[1] + psi[1]), log(2), rep(-20, 12), 0.5, 0.25)
    pars <- pm(x=x1, w=w1)

    capture.output(trees <- replicate(20, sim_bd_proc(n=40, l=pars$l, m=pars$m,
                                                      psi=pars$psi, init=1),
                                      simplify=FALSE))
    likwrap <- function(x, phylo) {
        w1[2] <- x[1]
        w1[15] <- x[2]
        w1[16] <- x[3]
        calc_bd_lm_nll(w=w1, x=x1, param_map=pm, y=list(phylo))
    }
    get_score <- function(phylo){
      numDeriv::grad(likwrap, x=w1[c(2, 15, 16)],
                     method="simple", phylo=phylo)
    }
    scores <- sapply(trees, get_score)
    htests <- apply(scores, 1, stats::t.test)
    ps <- sapply(htests, "[[", "p.value")
    expect_true(all(stats::p.adjust(ps) > 0.05))
})

test_that("Regularization path computed without error", {
    skip_on_cran()
    skip_if_not_installed("fizzlipuzzli")
    set.seed(2)

    pm <- gen_param_map(2, 1, .1)
    x1 <- cbind(c(1, 0, 0, 0))
    w1 <- c(0.73, 0.83, 0, 2)
    pars <- pm(x=x1, w=w1)

    capture.output(trees <- replicate(1, sim_bd_proc(n=10, l=pars$l, m=pars$m,
                                                      psi=pars$psi, init=1),
                                      simplify=FALSE))

    pf <- c(0, 0, rep(1, 2))
    init <- w1
    init[as.logical(pf)] <- 0
    out <- get_gpnet(x=x1, y=trees[1], calc_convex_nll=calc_bd_lm_nll,
                     param_map=pm, nlambda=100, lambda.min.ratio=0.5,
                     verbose=TRUE, penalty.factor=pf, thresh=1e-4,
                     winit=init, alpha=1)
    succeed()
})

test_that("gpnet estimates are reasonable", {
    skip_on_cran()
    skip_if_not_installed("fizzlipuzzli")
    set.seed(2)

    pm <- gen_param_map(2, 1, .1)
    x1 <- cbind(c(0, 1, 0, 0), c(0, 0, 1, 0), c(0, 0, 0, 1))
    w1 <- c(0.73, 0.83, 0, -4, 2, 1)
    pars <- pm(x=x1, w=w1)

    capture.output(trees <- replicate(1, sim_bd_proc(n=80, l=pars$l, m=pars$m,
                                                      psi=pars$psi, init=1),
                                      simplify=FALSE))

    pf <- c(0, 0, rep(1, 4))
    init <- w1
    init[as.logical(pf)] <- 0
    out <- get_gpnet(x=x1, y=trees[1], calc_convex_nll=calc_bd_lm_nll,
                     param_map=pm, nlambda=100, lambda.min.ratio=0.5,
                     verbose=TRUE, penalty.factor=pf, thresh=1e-4,
                     winit=init, alpha=1)
    succeed()
})

