set.seed(1)

context("phylogenetic likelihood")

test_that("Phylogenetic estimates are consistent", {
    skip_if_not_installed("phangorn")

    tree <- rcoal(4)
    data <- simSeq(tree, l = 1e5, type = "DNA", bf = c(.1, .2, .3, .4),
                   Q = rep(1, 6), rate=1)
    fitr <- pml(tree, data)

    fitr <- optim.pml(fitr, optBf = TRUE, optRooted = TRUE)
    stopifnot(all(tree$tip.label == fitr$tree$tip.label))
    ord1 <- order(tree$edge[, 1], tree$edge[, 2])
    ord2 <- order(fitr$tree$edge[, 1], fitr$tree$edge[, 2])

    expect_true(isTRUE(all.equal(tree$edge.length[ord1],
                                 fitr$tree$edge.length[ord2], tol = 0.01)))
    expect_true(isTRUE(all.equal(fitr$bf, 1:4 / 10, tol = 0.01)))
})

