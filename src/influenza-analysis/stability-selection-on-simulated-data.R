#!/usr/bin/Rscript

load("influenza-c3.RData")

sel_coef <- c(selected_fit$a0[, 10], selected_fit$beta[, 10])
sel_par <- pm1(xstable, sel_coef)
nsamples <- floor(length(tree_timel[[1]]$tip.label) / 2)

sim_tree <- penaltree::sim_bd_proc(n = nsamples, l = sel_par$l,
                                   m = sel_par$m, psi = sel_par$psi, init = 3)
while(max(penaltree:::get_nodeheights(sim_tree)$nh) > 8) {
  sim_tree <- penaltree::sim_bd_proc(n = nsamples, l = sel_par$l,
                                     m = sel_par$m, psi = sel_par$psi, init = 3)
}

init3 <- init1
sp_sim <- penaltree::stabpath_gpnet(x = x2, y = list(sim_tree),
               calc_convex_nll = penaltree::calc_bd_lm_nll,
               param_map = pm1, nlambda = 40, lambda.min.ratio = 0.1,
               make_log = TRUE, penalty.factor = pf1,
               thresh = 1e-3, winit = init3, alpha = 1,
               steps = 30, mc.cores = ncpu)
save.image("influenza-c4.RData")
