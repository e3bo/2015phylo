set.seed(1)

context("phylogenetic likelihood")

test_that("Phylogenetic estimates are consistent", {
    skip_if_not_installed("phangorn")
    skip_if_not_installed("ape")

    ntips <- 4
    tree <- ape::rcoal(ntips)
    data <- phangorn::simSeq(tree, l = 1e5, type = "DNA",
                             bf = c(.1, .2, .3, .4), Q = rep(1, 6), rate=1)
    fitr <- phangorn::pml(tree, data)

    fitr <- phangorn::optim.pml(fitr, optBf = TRUE, optRooted = TRUE)
    stopifnot(all(tree$tip.label == fitr$tree$tip.label))
    ord1 <- order(tree$edge[, 1], tree$edge[, 2])
    ord2 <- order(fitr$tree$edge[, 1], fitr$tree$edge[, 2])

    expect_true(isTRUE(all.equal(tree$edge.length[ord1],
                                 fitr$tree$edge.length[ord2], tol = 0.01)))
    expect_true(isTRUE(all.equal(fitr$bf, 1:4 / 10, tol = 0.01)))
})

test_that("Own optimization matches others", {
    skip_if_not_installed("phangorn")
    skip_if_not_installed("ape")

    ntips <- 10
    tree <- ape::rcoal(ntips)
    bf <- seq(1, 4) / 10
    data <- phangorn::simSeq(tree, l = 1e4, type = "DNA",
                             bf = bf, Q = rep(1, 6), rate=1)
    fitr <- phangorn::pml(tree, data, bf=bf)

    fitr <- phangorn::optim.pml(fitr, optRooted = TRUE)

    tiphts <- rep(0, ntips)

    obj <- function(x) {
        pml_wrapper(tree, nodeheights=x, tipheights=tiphts, data, bf=bf)
    }
    ans <- optim(seq(2, ntips)/10, obj, method="BFGS")

    ndtrue <- ape::node.depth.edgelength(tree)
    tree_est <- set_branchlengths(tree, nodeheights=ans$par, tipheights=tiphts)$tree
    ndest <- ape::node.depth.edgelength(tree_est)
    ndest2 <- ape::node.depth.edgelength(fitr$tree)

    expect_true(isTRUE(all.equal(ndest, ndest2, tol=.1)))
    expect_true(isTRUE(all.equal(-ans$val, as.numeric(logLik(fitr)), tol=1e-3)))
})

