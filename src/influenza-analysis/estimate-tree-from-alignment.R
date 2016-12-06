#!/usr/bin/Rscript

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
