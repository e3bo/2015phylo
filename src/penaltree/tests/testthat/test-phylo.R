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

test_that(paste("Able to estimate rate and rootheight with inferred topology"),{
    skip_if_not_installed("phangorn")
    skip_if_not_installed("ape")

    ntips <- 20
    tree <- ape::rtree(ntips)

    bf <- rep(.25, 4)
    rate <- 1e-3
    data <- phangorn::simSeq(tree, l = 1e3, type = "DNA",
                             bf = bf, Q = rep(1, 6), rate = rate)
    dists <- dist.ml(data)
    treeUPGMA <- upgma(dists)

    phyDat2msa <- function(data){
        alpha <- toupper(attr(data, "levels"))
        data <- toupper(as.character(data))
        cvl <- apply(data, 1, function(x) paste(x, collapse=""))
        rphast::msa(seqs = as.character(cvl), names = names(cvl),
                    alphabet = paste(alpha, collapse = ''))
    }
    msa <- phyDat2msa(data)
    treechar <- ape::write.tree(treeUPGMA)
    tmod <- rphast::tm(treechar, "JC69", backgd = rep(.25, 4))

    nh <- get_nodeheights(tree)
    obj <- function(x) {
        lmsa_wrapper(treeUPGMA, nodeheights = x[-1], tipheights = nh$tip, msa=msa, bf=bf, rate=x[1])
    }

    init <- c(0.1, runif(nh$node, max=3))
    init <- ifelse(init < 0, 0, init)
    ans <- rphast::optim.rphast(obj, init, lower=c(1e-6, rep(0, ntips-1)),
                                upper=c(1, rep(10, ntips -1)))#, logfile="/tmp/optim.log")
    tree_est <- set_branchlengths(treeUPGMA, nodeheights = ans$par[-1],
                                  tipheights = nh$tip)$tree
    nhest <- get_nodeheights(tree_est)$node

    expect_equal(sort(nhest), sort(nh$node), tol = .1)
    expect_equal(max(nhest), max(nh$node), tol=1)
    expect_equal(ans$par[1], rate, tol=1)
})

test_that(paste("Able to estimate parameters given HYK subs model",
                "with good starting point"),{
    skip_if_not_installed("phangorn")
    skip_if_not_installed("ape")

    ntips <- 20
    tree_time <- ape::rtree(ntips)

    kappa <- 4
    bf <- c(A=.25, C=.25, G=.1, T=.4)
    rate <- 1e-3

    Q <- get_hky_Q(kappa=4, pi=bf)
    subs_per_time <- 1e-3
    tree_subs <- tree_time
    tree_subs$edge.length <- tree_subs$edge.length * subs_per_time
    tree_char <- ape::write.tree(tree_subs)
    true_tm <- rphast::tm(tree_char, subst.mod="HKY85", rate.matrix=Q,
                          backgd=bf)
    ncols <- 1e3
    sim <- rphast::simulate.msa(true_tm, ncols)

    #likelihood.msa(sim, true_tm)
    nh <- get_nodeheights(tree_time)

    obj <- function(x) {
        subs_per_time <- x[1]
        pi <- x[seq(2, 5)]
        pi <- pi / sum(pi)
        names(pi) <- c("A", "C", "G", "T")
        kappa <- x[6]
        node_times <- x[seq(7, length(x))]
        lmsa_wrapper(tree_time, node_times = node_times, tip_times = nh$tip,
                     msa = sim, subs_per_time = subs_per_time,
                     subs_model = "HKY85",
                     subs_pars = list(kappa = kappa, pi = pi))
    }
    init <- c(subs_per_time, bf, 4, nh$node)
    init <- ifelse(init < 0, 0, init)
    ans <- rphast::optim.rphast(obj, init, lower=rep(0, length(init)),
                                logfile="/tmp/optim.log")
    nhest <- ans$par[-seq(1, 6)]
    tree_est <- set_branchlengths(tree_time, nodeheights = nhest,
                                  tipheights = nh$tip)$tree
    kappa_est <- ans$par[6]
    bf_est <- ans$par[seq(2, 5)]
    bf_est <- bf_est / sum(bf_est)

    expect_equal(nhest, nh$node, tol = .5)
    expect_equal(ans$par[1], subs_per_time, tol=.5)
    expect_equal(bf_est, unname(bf), tol=.5)
    expect_equal(kappa_est, kappa, tol=.5)
})

test_that(paste("Able to estimate parameters given HYK subs model",
                "with inferred topology"),{
    skip_if_not_installed("phangorn")
    skip_if_not_installed("ape")

    ntips <- 40
    tree_time <- ape::rtree(ntips)

    kappa <- 4
    bf <- c(A=.25, C=.25, G=.1, T=.4)
    rate <- 1e-3

    Q <- get_hky_Q(kappa=kappa, pi=bf)
    subs_per_time <- 1e-3
    tree_subs <- tree_time
    tree_subs$edge.length <- tree_subs$edge.length * subs_per_time

    tree_char <- ape::write.tree(tree_subs)
    true_tm <- rphast::tm(tree_char, subst.mod="HKY85", rate.matrix=Q,
                          backgd=bf)
    ncols <- 1e3
    sim <- rphast::simulate.msa(true_tm, ncols)
    charmat <- do.call(rbind, (strsplit(sim[[1]], split='')))
    rownames(charmat) <- sim$names
    simpd <- phyDat(charmat)
    write.phyDat(simpd, file="sim.phy")
    dist <- dist.ml(simpd, model="JC69")
    tree_upgma <- upgma(dist)

    #likelihood.msa(sim, true_tm)
    nh <- get_nodeheights(tree_time)

    obj <- function(x) {
        subs_per_time <- x[1]
        pi <- x[seq(2, 5)]
        pi <- pi / sum(pi)
        names(pi) <- c("A", "C", "G", "T")
        kappa <- x[6]
        node_times <- x[seq(7, length(x))]
        lmsa_wrapper(tree_upgma, node_times = node_times, tip_times = nh$tip,
                     msa = sim, subs_per_time = subs_per_time,
                     subs_model = "HKY85",
                     subs_pars = kappa, pi = pi)
    }
    init <- c(1e-4, rep(.25, 4), 2, runif(nh$node, max=3))
    init <- ifelse(init < 0, 0, init)
    ans <- rphast::optim.rphast(obj, init, lower=rep(0, length(init)),
                                 logfile="/tmp/optim.log")
    ## can repeat in case stopped due to max iterations, as noticeable only from looking at logfile
    #ans <- rphast::optim.rphast(obj, ans0$par, lower=rep(0, length(init)),
    #                            logfile="/tmp/optim.log")

                                        # optim allows iterations to be controlled but does much worse
    #ans2 <- optim(par=init, fn=obj, lower=rep(1e-6, length(init)),
    #              upper=c(1, rep(1, 4), 10, rep(6, 19)), method="L-BFGS-B", control=list(fnscale=-1, trace=5))

    nhest <- ans$par[-seq(1, 6)]
    tree_est <- set_branchlengths(tree_upgma, nodeheights = nhest,
                                  tipheights = nh$tip)$tree
    kappa_est <- ans$par[6]
    bf_est <- ans$par[seq(2, 5)]
    bf_est <- bf_est / sum(bf_est)

    expect_equal(sort(nhest), sort(nh$node), tol = .5)
    expect_equal(ans$par[1], subs_per_time, tol=.5)
    expect_equal(bf_est, unname(bf), tol=.5)
    expect_equal(kappa_est, kappa, tol=.5)
})

test_that(paste("Able to estimate parameters given HYK subs model + Gamma4",
                "with good starting point"),{
    skip_if_not_installed("phangorn")
    skip_if_not_installed("ape")

    ntips <- 50
    tree_time <- ape::rtree(ntips)

    kappa <- 4
    bf <- c(A=.25, C=.25, G=.1, T=.4)

    alpha <- 0.5
    nrates <- 4
    rate.consts <- phangorn::discrete.gamma(alpha, nrates)
    rate.weights <- rep(1 / nrates, nrates)

    subs_per_time <- 1e-1
    tree_subs <- tree_time
    tree_subs$edge.length <- tree_subs$edge.length * subs_per_time

    tree_char <- ape::write.tree(tree_subs)
    true_tm <- rphast::tm(tree_char, subst.mod="HKY85", backgd = bf,
                          nratecats = nrates, rate.consts = rate.consts,
                          rate.weights = rate.weights)
    true_tm <- rphast::set.rate.matrix.tm(true_tm, params=kappa)
    ncols <- 1e1
    sim <- rphast::simulate.msa(true_tm, ncols)

    nh <- get_nodeheights(tree_time)
    obj <- function(x) {
        subs_per_time <- x[1]
        pi <- x[seq(2, 5)]
        pi <- pi / sum(pi)
        names(pi) <- c("A", "C", "G", "T")
        kappa <- x[6]
        alpha <- x[7]
        node_times <- x[seq(8, length(x))]
        lmsa_wrapper(tree_time, node_times = node_times, tip_times = nh$tip,
                     msa = sim, subs_per_time = subs_per_time, alpha = alpha,
                     subs_model = "HKY85", nrates = 4,
                     subs_pars = kappa, pi = pi)
    }
    init <- c(subs_per_time, bf, kappa, alpha, nh$node)
    #init <- init * runif(init, max=3)
    init <- ifelse(init < 0, 0, init)
    ans <- rphast::optim.rphast(obj, init, lower=rep(0, length(init)),
                                 logfile="/tmp/optim.log")
    nhest <- ans$par[-seq(1, 7)]
    tree_est <- set_branchlengths(tree_time, nodeheights = nhest,
                                  tipheights = nh$tip)$tree
    kappa_est <- ans$par[6]
    alpha_est <- ans$par[7]
    bf_est <- ans$par[seq(2, 5)]
    bf_est <- bf_est / sum(bf_est)

    expect_equal(nhest, nh$node, tol = .5)
    expect_equal(ans$par[1], subs_per_time, tol=.5)
    expect_equal(bf_est, unname(bf), tol=.5)
    expect_equal(kappa_est, kappa, tol=.5)
    expect_equal(alpha_est, alpha, tol=.5)
})

test_that(paste("Able to estimate parameters given GTR subs model + Gamma4",
                "with inferred topology"),{
    skip_if_not_installed("phangorn")
    skip_if_not_installed("ape")

    ntips <- 50
    tree_time <- ape::rtree(ntips)
    nh <- get_nodeheights(tree_time)
    tree_time$tip.label <- paste(tree_time$tip.label,
                                 nh$tipheights[tree_time$tip.label], sep="_")
    names(nh$tipheights) <- paste(names(nh$tipheights), nh$tipheights, sep="_")

    subs_params <- c(AC=.2, AG=.8, AT=.14, CG=.9, CT=.1, GT=.12)
    subs_params <- subs_params / subs_params[6]
    bf <- c(.2, .14, .25, .4)
    bf <- bf / sum(bf)

    alpha <- 0.5
    nrates <- 4
    rate.consts <- phangorn::discrete.gamma(alpha, nrates)
    rate.weights <- rep(1 / nrates, nrates)

    subs_per_time <- 1e-1
    tree_subs <- tree_time
    tree_subs$edge.length <- tree_subs$edge.length * subs_per_time

    tree_char <- ape::write.tree(tree_subs)
    true_tm <- rphast::tm(tree_char, subst.mod="REV", backgd = bf,
                          nratecats = nrates, rate.consts = rate.consts,
                          rate.weights = rate.weights)
    true_tm <- rphast::set.rate.matrix.tm(true_tm, params=subs_params)
    ncols <- 1e3
    sim <- rphast::simulate.msa(true_tm, ncols)

    charmat <- do.call(rbind, (strsplit(sim[[1]], split='')))
    rownames(charmat) <- sim$names
    simpd <- phangorn::phyDat(charmat)
    simdnb <- ape::as.DNAbin(simpd)
    tr <- ips::raxml(simdnb, m="GTRGAMMA", p=12345, N=3, f="a",
                     exec="/usr/bin/raxmlHPC")
    bt <- tr$bestTree
    btr <- set_best_root(bt, nh$tipheights)
    rate_ests <- get_raxml_ests(tr=tr)
    temp_ests <- eval_temporal_signal(btr, nh$tipheights)

    obj <- function(x) {
        subs_per_time <- x[1]
        pi <- x[seq(2, 4)]
        pi <- c(pi, 1) / (1 + sum(pi))
        names(pi) <- c("A", "C", "G", "T")
        subs_pars <- c(x[seq(5, 9)], 1)
        alpha <- x[10]
        node_times <- x[seq(11, length(x))]
        lmsa_wrapper(btr, node_times = node_times, tip_times = nh$tip,
                     msa = sim, subs_per_time = subs_per_time, alpha = alpha,
                     subs_model = "REV", nrates = 4,
                     subs_pars = subs_pars, pi = pi)
    }
    init <- c(temp_ests$subs_per_time, rate_ests$bf[-4] / rate_ests$bf[4],
              rate_ests$gtr_pars[-6], rate_ests$alpha, nh$node)
    # TODO get initial node heights from temp ests
    ans <- rphast::optim.rphast(obj, init, lower=rep(0, length(init)),
                                 logfile="/tmp/optim.log")
    nhest <- ans$par[-seq(1, 10)]
    tree_est <- set_branchlengths(btr, nodeheights = nhest,
                                  tipheights = nh$tip)$tree
    subs_pars_est <-  c(ans$par[seq(5, 9)], 1)
    alpha_est <- ans$par[10]
    bf_est <- c(ans$par[seq(2, 4)], 1)
    bf_est <- bf_est / sum(bf_est)

    expect_lt(ape::dist.topo(tree_est, tree_time), length(tree_est$tip.label) - 3)
    expect_equal(sort(nhest), sort(nh$node), tol = .5)
    expect_equal(ans$par[1], subs_per_time, tol=.5)
    expect_equal(bf_est, unname(bf), tol=.5)
    expect_equal(subs_pars_est[-6], unname(subs_params)[-6], tol=.5)
    expect_equal(alpha_est, alpha, tol=.5)
})
