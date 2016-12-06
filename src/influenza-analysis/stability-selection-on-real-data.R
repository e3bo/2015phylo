#!/usr/bin/Rscript

load("influenza-c1.RData")

x2 <- diag(9)[, -1]
tree_timel <- lapply(tree_info, "[[", "tree_time")

pm1 <- penaltree::gen_param_map(3, ntrees=1, psampled=.1)
init1 <- c(1, -.10, 0, 0, rep(0, 8))
pf1 <- c(0, 0, rep(1, 10))

pars <- pm1(x = x2, w = init1)

tree_timel <- lapply(tree_info, "[[", "tree_time")
ncpu <- min(parallel::detectCores() - 1, 30)
sp <- penaltree::stabpath_gpnet(x = x2, y = tree_timel[1],
               calc_convex_nll = penaltree::calc_bd_lm_nll,
               param_map = pm1, nlambda = 10, lambda.min.ratio = 0.1,
               make_log = TRUE, penalty.factor = pf1,
               thresh = 1e-3, winit = init1, alpha = 1,
               steps=30, mc.cores = ncpu)
save.image("influenza-c2.RData")
