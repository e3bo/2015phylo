
library(ape)
raxmlbin <- "/usr/bin/raxmlHPC"
alignbin <- "../data/swine-influenza-alignments-dnabin.rds"

stopifnot(file.exists(raxmlbin))
stopifnot(file.exists(alignbin))

align <- readRDS(alignbin)

get_tr <- function(x) {
  ips::raxml(x, m = "GTRGAMMA", p = 12345, N = 3, f = "a", exec = raxmlbin)
}

tr <- lapply(align, get_tr)
rate_ests <- lapply(tr, get_raxml_ests)
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
    temp_ests <- eval_temporal_signal(btr, -td)
    nhinit <- get_time_tree_internal_nodeheights(btr, temp_ests$subs_per_time,
                                                 -td)

    metar <- strsplit(btr$tip.label, "_")
    btr$geo_states <- sapply(metar, "[[", 1)

    se_levs <- c("SC", "NC")
    mw_levs <- c("IL", "IN", "IA", "KS", "MI", "MO", "NE", "OH", "SD",
                 "WI")
    sc_levs <- c("TX", "OK", "MN")
    btr$states <- factor(btr$geo_states)
    levels(btr$states) <- list("se" = se_levs, "mw" = mw_levs, "sc" = sc_levs)
    btr$states <- as.integer(btr$states)

    tree_time <- set_branchlengths(btr, nhinit, -td)$tree
    list(td = td, temp_ests = temp_ests, btr = btr, nhinit = nhinit,
         tree_time = tree_time)    
}
tree_info <- lapply(bt, process_trees)

pm <- gen_param_map(3, ntrees=length(tree_info))
init <- c(0, 0, 0, 0, 0, 0, rep(0, 8))

x2 <- diag(9)[, -1]

pars <- pm(x=x2, w=init)

tree_timel <- lapply(tree_info, "[[", "tree_time")
out <- get_gpnet(x=x2, y=tree_timel, calc_convex_nll=calc_bd_lm_nll,
                 param_map=pm, nlambda=13, lambda.min.ratio=0.75,
                 verbose=TRUE, penalty.factor=c(0, 0, rep(1,12)),
                 thresh=1e-4, winit=init, alpha=1)

pm1 <- gen_param_map(3, ntrees=1, psampled=.1)
init1 <- c(0, 0, 0, 0, rep(0, 8))
pf1 <- c(0, 0, rep(1, 10))

pars <- pm1(x = x2, w = init1)

tree_timel <- lapply(tree_info, "[[", "tree_time")
out1 <- get_gpnet(x = x2, y = tree_timel[1], calc_convex_nll = calc_bd_lm_nll,
                  param_map = pm1, nlambda = 13, lambda.min.ratio = 0.75,
                  verbose = TRUE, penalty.factor = pf1,
                  thresh = 1e-4, winit = init1, alpha = 1)

out2 <- get_gpnet(x = x2, y = tree_timel[2], calc_convex_nll = calc_bd_lm_nll,
                  param_map = pm1, nlambda = 13, lambda.min.ratio = 0.75,
                  verbose = TRUE, penalty.factor = pf1,
                  thresh = 1e-4, winit = init1, alpha = 1)

obj <- function(x) {
    foo[c(1,2)] <- x
    calc_bd_lm_nll(w=foo, x=x2, y=tree_timel[1], param_map=pm1)
}
ans <- optim(c(-1, -1), obj)





calc_bd_lm_nll(foo, x2, tree_timel[1], pm1)
