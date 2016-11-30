#!/usr/bin/Rscript

raxmlbin <- "/usr/bin/raxmlHPC"
alignbin <- "../data/swine-influenza-alignments-dnabin.rds"

stopifnot(file.exists(raxmlbin))
stopifnot(file.exists(alignbin))

align <- readRDS(alignbin)

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

pm <- penaltree::gen_param_map(3, ntrees=length(tree_info))
init <- c(0, 0, 0, 0, 0, 0, rep(0, 8))

x2 <- diag(9)[, -1]

pars <- pm(x=x2, w=init)

tree_timel <- lapply(tree_info, "[[", "tree_time")

pm1 <- penaltree::gen_param_map(3, ntrees=1, psampled=.1)
init1 <- c(1, -.10, 0, 0, rep(0, 8))
pf1 <- c(0, 0, rep(1, 10))

pars <- pm1(x = x2, w = init1)

tree_timel <- lapply(tree_info, "[[", "tree_time")
out1 <- penaltree::get_gpnet(x = x2, y = tree_timel[1],
                  calc_convex_nll = penaltree::calc_bd_lm_nll,
                  param_map = pm1, nlambda = 50, lambda.min.ratio = 0.5,
                  verbose = TRUE, penalty.factor = pf1,
                  thresh = 1e-4, winit = init1, alpha = 1)

save.image("influenza.RData")

q('no')

out <- get_gpnet(x = x2, y = tree_timel, calc_convex_nll=penaltree::calc_bd_lm_nll,
                 param_map=pm, nlambda=13, lambda.min.ratio=0.75,
                 verbose=TRUE, penalty.factor=c(0, 0, rep(1,12)),
                 thresh=1e-4, winit=init, alpha=1)

out2 <- penaltree::get_gpnet(x = x2, y = tree_timel[2], calc_convex_nll = calc_bd_lm_nll,
                  param_map = pm1, nlambda = 13, lambda.min.ratio = 0.75,
                  verbose = TRUE, penalty.factor = pf1,
                  thresh = 1e-4, winit = init1, alpha = 1)

obj <- function(x) {
    foo[c(1,2)] <- x
    l <- calc_bd_lm_nll(w=foo, x=x2, y=tree_timel[1], param_map=pm1)
    print(paste(c(x, l)))
    l
}
ans <- optim.rphast(obj, c(1, -.1), lower=c(-1,-1), upper=c(1, 0.5), logfile="/tmp/a.log")



p <- tree_timel[[2]]
p <- addroot(p, 0)
p$states[p$states == 3] <- 2
 par <- c(1, 1, 1, 1)
 fix <- rbind(c(1, 6, 7, 8), c(1, -5, 0, 0), c(1, 1, 1, 1))
 tplik <- TreePar::LikTypesSTT(par=par, phylo=p, fix=fix, sampfrac=s,
                                  survival=0, freq=0.1, posR=0)



calc_bd_lm_nll(foo, x2, tree_timel[1], pm1)
