
raxmlbin <- "/usr/bin/raxmlHPC"
stopifnot(file.exists(raxmlbin))

nsi <- ape::read.dna("../work-sshd/nonsIndel-aligned.fasta-gb",
                     format = "fasta")

OH_variant_ind <- grep("OH_14-Jun", dimnames(nsi)[[1]])
nsi2 <- nsi[-OH_variant_ind, ]

tr <- ips::raxml(nsi2, m = "GTRGAMMA", p = 12345, N = 3, f = "a",
                 exec = raxmlbin)
bt <- tr$bestTree

meta <- strsplit(bt$tip.label, "_")
tip_dates <- as.Date(sapply(meta, "[[", 2), format = "%d-%b-%Y")
tip_states <- sapply(meta, "[[", 1)
td <- as.integer(tip_dates - max(tip_dates)) / 365
names(td) <- bt$tip.label

btr <- ape::rtt(bt, tip.dates = td, objective = "rms")
temp_ests <- penaltree::eval_temporal_signal(btr, -td)

metar <- strsplit(btr$tip.label, "_")
btr$geo_states <- sapply(metar, "[[", 1)
btr$states <- as.integer(btr$geo_states %in% c("IA", "MN")) + 1

# todrop <- sample.int(length(btr$tip.label), 1)
# btrs <- drop.tip(btr, todrop)
# btrs$states <- btrs$states[-todrop]
#
# nhinit <- get_time_tree_internal_nodeheights(btrs,
#                                  temp_ests$subs_per_time, -td[-todrop])
# tree_time <- set_branchlengths(btrs, nhinit, -td[-todrop])$tree

nhinit <- penaltree::get_time_tree_internal_nodeheights(btr,
                                                        temp_ests$subs_per_time,
                                                        -td)
tree_time <- penaltree::set_branchlengths(btr, nhinit, -td)$tree

pm <- penaltree::gen_param_map(2, 1, .1)
init <- c(0.73, 0.83, 0, 0, 0, 0)
x2 <- cbind(c(0, 0, 0, 1), c(1, 0, 0, 0), c(0, 1, 0, 0))
pf <- c(0, 0, 1, 1, 1, 1)

pars <- pm(x = x2, w = init)
out <- penaltree::get_gpnet(x = x2, y = list(tree_time),
                            calc_convex_nll = calc_bd_lm_nll, param_map = pm,
                            nlambda = 100, lambda.min.ratio = 0.5,
                            verbose = TRUE, penalty.factor = pf, thresh = 1e-4,
                            winit = init, alpha = 1)
save.image("pedv-wg.RData")
