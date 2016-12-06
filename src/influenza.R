#!/usr/bin/Rscript

library(c060)
set.seed(1)
raxmlbin <- "/usr/bin/raxmlHPC"
stopifnot(file.exists(raxmlbin))

align <- list(H1N1 = penaltree::H1N1_alignment)

get_tr <- function(x) {
  ips::raxml(x, m = "GTRGAMMA", p = 12345, N = 3, f = "a", exec = raxmlbin)
}

tr <- lapply(align, get_tr)
rate_ests <- lapply(tr, penaltree:::get_raxml_ests)
bt <- lapply(tr, "[[", "bestTree")

process_trees <- function(tree){
    meta <- strsplit(tree$tip.label, "_")
    tmpf <- function(x){
        x[[length(x)]]
    }
    td <- as.numeric(sapply(meta, tmpf))
    td <- td - max(td)
    names(td) <- tree$tip.label
    btr <- ape::rtt(tree, tip.dates = td, objective = "rms")
    temp_ests <- penaltree::eval_temporal_signal(btr, -td)
    nhinit <- penaltree::get_time_tree_internal_nodeheights(btr,
                                                  temp_ests$subs_per_time, -td)

    metar <- strsplit(btr$tip.label, "_")
    btr$geo_states <- sapply(metar, "[[", 1)

    se_levs <- c("SC", "NC")
    mw_levs <- c("MN", "WI", "IA")
    other_levs <- c("IL", "IN", "KS", "MI", "MO", "NE", "OH", "SD")
    btr$states <- factor(btr$geo_states)
    levels(btr$states) <- list("se" = se_levs, "mw" = mw_levs, "other" = other_levs)
    btr$states <- as.integer(btr$states)

    tree_time <- penaltree::set_branchlengths(btr, nhinit, -td)$tree
    list(td = td, temp_ests = temp_ests, btr = btr, nhinit = nhinit,
         tree_time = tree_time)
}
tree_info <- lapply(bt, process_trees)
save.image("influenza-c1.RData")

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

spstats <- plot(sp)
xstable <- x2[, 4L - 2L, drop=FALSE]
stopifnot("V4" %in% names(spstats$stable))

pf2 <- c(0, 0, 100, 100, 1)
init2 <- init1[c(seq(1, 4), 4L + 4L)]
selected_fit <- penaltree::get_gpnet(x = xstable, y = tree_timel[1],
               calc_convex_nll = penaltree::calc_bd_lm_nll,
               param_map = pm1, nlambda = 10, lambda.min.ratio = 0.01,
               make_log = TRUE, penalty.factor = pf2,
               thresh = 1e-3, winit = init2, alpha = 1)
save.image("influenza-c3.RData")

sel_coef <- c(selected_fit$a0[, 10], selected_fit$beta[, 10])
sel_par <- pm1(xstable, sel_coef)
nsamples <- floor(length(tree_timel[[1]]$tip.label) * 2 / 3)

sim_tree <- penaltree::sim_bd_proc(n = nsamples, l = sel_par$l,
                                   m = sel_par$m, psi = sel_par$psi, init = 3)

# do stability selection and fitting of simulated tree

init3 <- init1
sp_sim <- penaltree::stabpath_gpnet(x = x2, y = list(sim_tree),
               calc_convex_nll = penaltree::calc_bd_lm_nll,
               param_map = pm1, nlambda = 20, lambda.min.ratio = 0.1,
               make_log = TRUE, penalty.factor = pf1,
               thresh = 1e-3, winit = init3, alpha = 1,
               steps = 30, mc.cores = ncpu)
save.image("influenza-c4.RData")
spstats_sim <- plot(sp_sim)

stopifnot(all(as.integer(spstats_sim$stable) == 4L)))

sim_fit <- penaltree::get_gpnet(x = xstable, y = list(sim_tree),
               calc_convex_nll = penaltree::calc_bd_lm_nll,
               param_map = pm1, nlambda = 10, lambda.min.ratio = 0.01,
               make_log = TRUE, penalty.factor = pf2,
               thresh = 1e-3, winit = init2, alpha = 1)

save.image("influenza-c5.RData")
q("no")
