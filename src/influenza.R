
library(ape)
raxmlbin <- "/usr/bin/raxmlHPC"
alignbin <- "../data/swine-influenza-alignments-dnabin.rds"

stopifnot(file.exists(raxmlbin))
stopifnot(file.exists(alignbin))


align <- readRDS(alignbin)

tr <- ips::raxml(align$H1N1, m = "GTRGAMMA", p = 12345, N = 3, f = "a", exec = raxmlbin)
rate_ests <- get_raxml_ests(tr = tr)
bt <- tr$bestTree

meta <- strsplit(bt$tip.label, "_")
td <- as.numeric(sapply(meta, "[[", 3))
td <- td - max(td)
names(td) <- bt$tip.label

btr <- ape::rtt(bt, tip.dates = td, objective = "rms")
temp_ests <- eval_temporal_signal(btr, -td)

metar <- strsplit(btr$tip.label, "_")
btr$geo_states <- sapply(metar, "[[", 1)
btr$states <- ifelse(btr$geo_states == "NC", 1, 2)
btr$states <- ifelse(btr$geo_states == "MN", 3, btr$states)

nhinit <- get_time_tree_internal_nodeheights(btr, temp_ests$subs_per_time, -td)
tree_time <- set_branchlengths(btr, nhinit, -td)$tree

pm <- gen_param_map(3)
init <- c(0, 0, 0.33, 0.33, rep(0, 8))

x2 <- diag(9)[, -1]

pars <- pm(x=x2, w=init)
out <- get_gpnet(x=x2, y=tree_time, calc_convex_nll=calc_bd_lm_nll,
                 param_map=pm, nlambda=100, lambda.min.ratio=0.1,
                 verbose=TRUE, penalty.factor=c(0, 0, rep(1,10)),
                 thresh=1e-3, winit=init, alpha=1)
