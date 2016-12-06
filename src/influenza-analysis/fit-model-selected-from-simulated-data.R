#!/usr/bin/Rscript

library(c060)
load("influenza-c4.RData")
spstats_sim <- plot(sp_sim)

stopifnot(all(as.integer(spstats_sim$stable) == 4L))

sim_fit <- penaltree::get_gpnet(x = xstable, y = list(sim_tree),
               calc_convex_nll = penaltree::calc_bd_lm_nll,
               param_map = pm1, nlambda = 10, lambda.min.ratio = 0.01,
               make_log = TRUE, penalty.factor = pf2,
               thresh = 1e-3, winit = init2, alpha = 1)

save.image("influenza-c5.RData")
