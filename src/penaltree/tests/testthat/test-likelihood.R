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
    ptlik <- calc_bdlik(l=l, m=m, psi=psi, freq=0.1, phylo=tree, survival=FALSE)
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
    ptlik <- calc_bdlik(l=l, m=m, psi=psi, freq=c(0.9, 0.1), phylo=tree,
                        survival=FALSE)
    expect_equal(-36.2374678051509, ptlik)

    ntypes <- 14
    l <- matrix(2, ncol=ntypes, nrow=ntypes)
    m <- rep(1, ntypes)
    psi <- rep(0.05, ntypes)
    ptlik <- calc_bdlik(l=l, m=m, psi=psi, freq=rep(1/ntypes, ntypes - 1),
                        phylo=tree, survival=FALSE)
    expect_equal(40.4926005495144, ptlik)
})
