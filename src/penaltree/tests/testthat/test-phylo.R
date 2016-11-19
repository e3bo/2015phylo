set.seed(1)

context("phylogenetic likelihood")

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

test_that(paste("Able to estimate parameters given GTR subs model + Gamma4",
                "with inferred topology"),{
    skip_if_not_installed("phangorn")
    skip_if_not_installed("ape")
    skip_if_not_installed("ips")
    raxmlbin <- "/usr/bin/raxmlHPC"
    skip_if_not(file.exists(raxmlbin))
    skip_on_cran()

    load("testdata.rda")

    x <- x[seq(1, 169), ]
    x1 <- x[, c(1, 2), drop=FALSE]

    pm <- gen_param_map(13)
    w1 <- c(log(2), 0.5, 0.25)
    pars <- pm(x=x1, w=w1)

    ntips <- 40
    capture.output(trees <- replicate(1, sim_bd_proc(n=ntips, l=pars$l, m=pars$m,
                                                      psi=pars$psi, init=1),
                                      simplify=FALSE))
    tree_time <- trees[[1]]
    nh <- get_nodeheights(tree_time)
    tree_time$tip.label <- paste(tree_time$tip.label,
                                 nh$tipheights[tree_time$tip.label], sep = "_")
    names(nh$tipheights) <- paste(names(nh$tipheights), nh$tipheights,
                                  sep = "_")

    subs_params <- c(AC = .2, AG = .8, AT = .14, CG = .9, CT = .1, GT = .12)
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
    true_tm <- rphast::tm(tree_char, subst.mod = "REV", backgd = bf,
                          nratecats = nrates, rate.consts = rate.consts,
                          rate.weights = rate.weights)
    true_tm <- rphast::set.rate.matrix.tm(true_tm, params = subs_params)
    ncols <- 1e3
    sim <- rphast::simulate.msa(true_tm, ncols)

    charmat <- do.call(rbind, (strsplit(sim[[1]], split = '')))
    rownames(charmat) <- sim$names
    simpd <- phangorn::phyDat(charmat)
    simdnb <- ape::as.DNAbin(simpd)
    tr <- ips::raxml(simdnb, m = "GTRGAMMA", p = 12345, N = 3, f = "a",
                     exec = raxmlbin)
    bt <- tr$bestTree
    tip.dates <- -nh$tipheights[bt$tip.label]
    btr <- ape::rtt(bt, tip.dates=tip.dates, objective="rms")

    rate_ests <- get_raxml_ests(tr = tr)
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
    nhinit <- get_time_tree_internal_nodeheights(btr, temp_ests$subs_per_time,
                                                 nh$tip)
    init <- c(temp_ests$subs_per_time, rate_ests$bf[-4] / rate_ests$bf[4],
              rate_ests$gtr_pars[-6], rate_ests$alpha, nhinit)
    (logfile <- tempfile(fileext = ".log"))
    ans <- rphast::optim.rphast(obj, init, lower = rep(0, length(init)),
                                 logfile = logfile)
    nhest <- ans$par[-seq(1, 10)]
    tree_est <- set_branchlengths(btr, nodeheights = nhest,
                                  tipheights = nh$tip)$tree
    subs_pars_est <-  c(ans$par[seq(5, 9)], 1)
    alpha_est <- ans$par[10]
    bf_est <- c(ans$par[seq(2, 4)], 1)
    bf_est <- bf_est / sum(bf_est)

    expect_lt(ape::dist.topo(tree_est, tree_time),
              length(tree_est$tip.label) - 3)
    expect_equal(sort(nhest), sort(nh$node), tol = .5)
    expect_equal(ans$par[1], subs_per_time, tol = .5)
    expect_equal(bf_est, unname(bf), tol = .5)
    expect_equal(subs_pars_est[-6], unname(subs_params)[-6], tol = .5)
    expect_equal(alpha_est, alpha, tol = .5)
})

test_that(paste("Able to estimate parameters given GTR subs model + Gamma4",
                "with inferred topology using own optimizer"),{
    skip_if_not_installed("phangorn")
    skip_if_not_installed("ape")
    skip_if_not_installed("ips")
    raxmlbin <- "/usr/bin/raxmlHPC"
    skip_if_not(file.exists(raxmlbin))
    skip_on_cran()

    load("testdata.rda")

    x <- x[seq(1, 169), ]
    x1 <- x[, c(1, 2), drop=FALSE]

    pm <- gen_param_map(13)
    w1 <- c(log(2), 0.5, 0.25)
    pars <- pm(x=x1, w=w1)

    ntips <- 40
    capture.output(trees <- replicate(1, sim_bd_proc(n=ntips, l=pars$l, m=pars$m,
                                                      psi=pars$psi, init=1),
                                      simplify=FALSE))
    tree_time <- trees[[1]]
    nh <- get_nodeheights(tree_time)
    tree_time$tip.label <- paste(tree_time$tip.label,
                                 nh$tipheights[tree_time$tip.label], sep = "_")
    names(nh$tipheights) <- paste(names(nh$tipheights), nh$tipheights,
                                  sep = "_")

    subs_params <- c(AC = .2, AG = .8, AT = .14, CG = .9, CT = .1, GT = .12)
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
    true_tm <- rphast::tm(tree_char, subst.mod = "REV", backgd = bf,
                          nratecats = nrates, rate.consts = rate.consts,
                          rate.weights = rate.weights)
    true_tm <- rphast::set.rate.matrix.tm(true_tm, params = subs_params)
    ncols <- 1e3
    sim <- rphast::simulate.msa(true_tm, ncols)

    charmat <- do.call(rbind, (strsplit(sim[[1]], split = '')))
    rownames(charmat) <- sim$names
    simpd <- phangorn::phyDat(charmat)
    simdnb <- ape::as.DNAbin(simpd)
    tr <- ips::raxml(simdnb, m = "GTRGAMMA", p = 12345, N = 3, f = "a",
                     exec = raxmlbin)
    bt <- tr$bestTree
    tip.dates <- -nh$tipheights[bt$tip.label]
    btr <- ape::rtt(bt, tip.dates = tip.dates, objective = "rms")

    rate_ests <- get_raxml_ests(tr = tr)
    temp_ests <- eval_temporal_signal(btr, nh$tipheights)

    pm <- gen_param_map_phylo_bd(btr, nh$tip)
    nhinit <- get_time_tree_internal_nodeheights(btr, temp_ests$subs_per_time,
                                                 nh$tip)
    init <- c(temp_ests$subs_per_time, rate_ests$bf[-4] / rate_ests$bf[4],
              rate_ests$gtr_pars[-6], rate_ests$alpha, nhinit, 0)
    (logfile <- tempfile(fileext = ".log"))

    pf <- c(rep(0, length(init) - 1), 1)
    out <- get_gpnet(x = matrix(1, ncol = 1), y = sim,
                     calc_convex_nll = calc_phylo_nll_bd, param_map = pm,
                     nlambda = 100, lambda.min.ratio = 0.1, verbose = TRUE,
                     winit = init, penalty.factor = pf, thresh = 1e-3)
    nhest <- unname(out$a0[-seq(1, 10), 1])
    tree_est <- set_branchlengths(btr, nodeheights = nhest,
                                  tipheights = nh$tip)$tree
    subs_pars_est <-  unname(c(out$a0[seq(5, 9), 1], 1))
    alpha_est <- unname(out$a0[10, 1])
    bf_est <- c(out$a0[seq(2, 4), 1], 1)
    bf_est <- unname(bf_est / sum(bf_est))
    est_subs_per_time <- unname(out$a0[1,1])

    expect_lt(ape::dist.topo(tree_est, tree_time),
              length(tree_est$tip.label) - 3)
    expect_equal(sort(nhest), sort(nh$node), tol = .5)
    expect_equal(est_subs_per_time, subs_per_time, tol = .5)
    expect_equal(bf_est, unname(bf), tol = .5)
    expect_equal(subs_pars_est[-6], unname(subs_params)[-6], tol = .5)
    expect_equal(alpha_est, alpha, tol = .5)
})
