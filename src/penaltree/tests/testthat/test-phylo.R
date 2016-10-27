set.seed(1)

context("phylogenetic likelihood")

test_that("Phylogenetic estimates are consistent", {
    skip_if_not_installed("phangorn")
    skip_if_not_installed("ape")

    tree <- ape::rcoal(4)
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

    tree <- ape::rcoal(4)
    bf <- seq(1, 4) / 10
    data <- phangorn::simSeq(tree, l = 1e5, type = "DNA",
                             bf = bf, Q = rep(1, 6), rate=1)
    obj <- function(x) {
        pml_wrapper(tree, nodeheights=x, tipheights=rep(0, 4), data, bf=bf)
    }
    ans <- optim(c(.2, .3, .4), obj, method="BFGS")

    ndtrue <- ape::node.depth.edgelength(tree)
    tree_est <- set_branchlengths(tree, nodeheights=ans$par, tipheights=rep(0,4))$tree
    ndest <- ape::node.depth.edgelength(tree_est)

    all.equal(ndest, ndtrue)

    par(mfrow=c(2,1))
    plot(tree)
    plot(tree_est)


})

