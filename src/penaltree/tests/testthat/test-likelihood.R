set.seed(3)

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
