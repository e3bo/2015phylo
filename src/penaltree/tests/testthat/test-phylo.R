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

    ntips <- 40
    tree <- ape::rcoal(ntips)
    bf <- rep(.25, 4)
    data <- phangorn::simSeq(tree, l = 1000, type = "DNA",
                             bf = bf, Q = rep(1, 6), rate = 1)
    fitr <- phangorn::pml(tree, data, bf = bf)

    fitr <- phangorn::optim.pml(fitr, optRooted = TRUE)

    tiphts <- rep(0, ntips)

    ndtrue <- ape::node.depth.edgelength(tree)
    ndest2 <- ape::node.depth.edgelength(fitr$tree)

    phyDat2msa <- function(data){
        alpha <- toupper(attr(data, "levels"))
        data <- toupper(as.character(data))
        cvl <- apply(data, 1, function(x) paste(x, collapse=""))
        rphast::msa(seqs = as.character(cvl), names = names(cvl),
                    alphabet = paste(alpha, collapse = ''))
    }
    msa <- phyDat2msa(data)
    treechar <- ape::write.tree(tree)
    tmod <- rphast::tm(treechar, "JC69", backgd = rep(.25, 4))
    rphast::likelihood.msa(x = msa, tm = tmod)
    pf <- rphast::phyloFit(msa = msa, init.mod = tmod, clock = TRUE)
    ptree <- ape::read.tree(text = pf$tree)
    ndest3 <- ape::node.depth.edgelength(ptree)

    expect_true(isTRUE(all.equal(ndest2, ndest3, tol = .1)))
    expect_true(isTRUE(all.equal(pf$likelihood, as.numeric(logLik(fitr)),
                                 tol = 1e-3)))

    obj <- function(x) {
        pml_wrapper(tree, nodeheights = x, tipheights = tiphts, data, bf = bf)
    }

    tnh <- get_nodeheights(tree)
    ans <- optim(tnh$nodeheights, obj, method = "L-BFGS-B")
    tree_est <- set_branchlengths(tree, nodeheights = ans$par,
                                  tipheights = tiphts)$tree
    ndest <- ape::node.depth.edgelength(tree_est)

    expect_true(isTRUE(all.equal(ndest, ndest2, tol = .1)))
    expect_true(isTRUE(all.equal(-ans$val, as.numeric(logLik(fitr)),
                                 tol = 1e-3)))
})
