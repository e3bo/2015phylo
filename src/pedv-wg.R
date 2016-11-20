
library(ape)
raxmlbin <- "/usr/bin/raxmlHPC"
stopifnot(file.exists(raxmlbin))

nsi <- read.dna("../work-sshd/nonsIndel-aligned.fasta-gb", format="fasta")

OH_variant_ind <- grep("OH_14-Jun", dimnames(nsi)[[1]])
nsi2 <- nsi[-OH_variant_ind, ]

tr <- ips::raxml(nsi2, m = "GTRGAMMA", p = 12345, N = 3, f = "a", exec = raxmlbin)
rate_ests <- get_raxml_ests(tr = tr)
bt <- tr$bestTree

meta <- strsplit(bt$tip.label, "_")
tip_dates <- as.Date(sapply(meta, "[[", 2), format="%d-%b-%Y")
tip_states <- sapply(meta, "[[", 1)
td <- as.integer(tip_dates - max(tip_dates))
names(td) <- bt$tip.label

# if we want to check with TempEst
#write.tree(bt, file="bt.tree")

btr <- ape::rtt(bt, tip.dates = td, objective = "rms")
temp_ests <- eval_temporal_signal(btr, -td)

metar <- strsplit(btr$tip.label, "_")
btr$states <- sapply(metar, "[[", 1)

nhinit <- get_time_tree_internal_nodeheights(btr, temp_ests$subs_per_time, -td)
tree_time <- set_branchlengths(btr, nhinit, -td)$tree
