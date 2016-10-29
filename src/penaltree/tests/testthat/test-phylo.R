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

    ntips <- 20
    tree <- ape::rcoal(ntips)
    tiphts <- rep(0, ntips)
    bf <- rep(.25, 4)
    data <- phangorn::simSeq(tree, l = 1000, type = "DNA",
                             bf = bf, Q = rep(1, 6), rate = 1)

    fitr <- phangorn::pml(tree, data, bf = bf)
    fitr <- phangorn::optim.pml(fitr, optRooted = TRUE)

    nhest2 <- get_nodeheights(fitr$tree)$node

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
    nhest3 <- get_nodeheights(ptree)$node

    expect_true(isTRUE(all.equal(nhest2, nhest3, tol = .1, scale=1)))
    expect_true(isTRUE(all.equal(pf$likelihood, as.numeric(logLik(fitr)),
                                 tol = 1e-3)))

    obj <- function(x) {
        lmsa_wrapper(tree, nodeheights = x, tipheights = tiphts, msa=msa,
                     bf=bf)
    }

    tnh <- get_nodeheights(tree)$node
    ans <- rphast::optim.rphast(obj, tnh, lower=rep(0, ntips - 1),
                                upper=rep(6, ntips - 1))
    tree_est <- set_branchlengths(tree, nodeheights = ans$par,
                                  tipheights = tiphts)$tree
    nhest <- get_nodeheights(tree_est)$node

    #plot(data.frame(tnh, nhest, nhest2, nhest3))
    expect_true(isTRUE(all.equal(nhest, nhest3, tol = .04, scale=1)))
    expect_true(isTRUE(all.equal(-ans$val, as.numeric(logLik(fitr)),
                                 tol = 1e-3)))
})

test_that("Node height calculation is correct", {

    tree <- ape::rtree(30)
    nde <- ape::node.depth.edgelength(tree)
    nh <- get_nodeheights(tree)
    ndc <- max(nh$nh) - nh$nh
    expect_equal(sort(ndc), sort(nde))

    tree2 <- set_branchlengths(tree, nh$nodeheights, nh$tipheights)$tree
    nde2 <- ape::node.depth.edgelength(tree2)
    expect_equal(nde, nde2)
})

test_that("With serial sampling, own optimization finds optimum near truth", {
    skip_if_not_installed("phangorn")
    skip_if_not_installed("ape")

    ntips <- 10
    tree <- ape::rtree(ntips)

    bf <- rep(.25, 4)
    data <- phangorn::simSeq(tree, l = 1000, type = "DNA",
                             bf = bf, Q = rep(1, 6), rate = 1)

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

    nh <- get_nodeheights(tree)
    obj <- function(x) {
        lmsa_wrapper(tree, nodeheights = x, tipheights = nh$tip, msa=msa, bf=bf)
    }

    ans <- rphast::optim.rphast(obj, nh$node, lower=rep(0, ntips - 1),
                                upper=rep(6, ntips - 1))
    tree_est <- set_branchlengths(tree, nodeheights = ans$par,
                                  tipheights = nh$tip)$tree
    nhest <- get_nodeheights(tree_est)$node

    expect_true(isTRUE(all.equal(nhest, nh$node, tol = .05, scale=1)))
    expect_true(isTRUE(all.equal(-ans$val, as.numeric(logLik(fitr)),
                                 tol = 1e-3)))
})

test_that(paste("With serial sampling, able to find decent estimate",
                "when given decent starting point"), {
    skip_if_not_installed("phangorn")
    skip_if_not_installed("ape")

    ntips <- 20
    tree <- ape::rtree(ntips)

    bf <- rep(.25, 4)
    data <- phangorn::simSeq(tree, l = 1000, type = "DNA",
                             bf = bf, Q = rep(1, 6), rate = 1)

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

    nh <- get_nodeheights(tree)
    obj <- function(x) {
        lmsa_wrapper(tree, nodeheights = x, tipheights = nh$tip, msa=msa, bf=bf)
    }

    init <- nh$node + runif(nh$node, min=-1, max=1)
    init <- ifelse(init < 0, 0, init)
    ans <- rphast::optim.rphast(obj, init, lower=rep(0, ntips - 1), upper=rep(6, ntips - 1))
    tree_est <- set_branchlengths(tree, nodeheights = ans$par,
                                  tipheights = nh$tip)$tree
    nhest <- get_nodeheights(tree_est)$node

    expect_true(isTRUE(all.equal(nhest, nh$node, tol = .05, scale=1)))
})

test_that(paste("Able to jointly estimate branchlengths and rate"), {
    skip_if_not_installed("phangorn")
    skip_if_not_installed("ape")

    ntips <- 20
    tree <- ape::rtree(ntips)

    bf <- rep(.25, 4)
    rate <- 1e-3
    data <- phangorn::simSeq(tree, l = 1e3, type = "DNA",
                             bf = bf, Q = rep(1, 6), rate = rate)

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

    nh <- get_nodeheights(tree)
    obj <- function(x) {
        lmsa_wrapper(tree, nodeheights = x[-1], tipheights = nh$tip, msa=msa, bf=bf, rate=x[1])
        # lmsa_wrapper(tree, nodeheights = x, tipheights = nh$tip, msa=msa, bf=bf, rate=1)
    }

    init <- c(0.1, nh$node + runif(nh$node, min=-1, max=1))
    init <- ifelse(init < 0, 0, init)
    ans <- rphast::optim.rphast(obj, init, lower=c(1e-6, rep(0, ntips-1)),
                                upper=c(1, rep(6, ntips -1)))#, logfile="/tmp/optim.log")
    tree_est <- set_branchlengths(tree, nodeheights = ans$par[-1],
                                  tipheights = nh$tip)$tree
    nhest <- get_nodeheights(tree_est)$node

    expect_true(isTRUE(all.equal(nhest, nh$node, tol = .5, scale=1)))
    expect_equal(ans$par[1], rate, tol=1)
})
