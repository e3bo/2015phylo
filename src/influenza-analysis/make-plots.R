load("influenza-c5.RData")

library(maps)
map(database = "state")

levs <- list(se=c("SC", "NC"), mw=c("MN", "WI", "IA"),
             ot=c("IL", "IN", "KS", "MI", "MO", "NE", "OH", "SD"))

abb2name <- function(abb){
    state.name[match(abb, state.abb)]
}
regions <- lapply(levs, abb2name)

png("region-map.png", width=800, height=800)
map(database = "state")
map(database = "state",regions = regions$se , col = "blue", fill=T, add=TRUE)
map(database = "state",regions = regions$mw , col = "orange", fill=T, add=TRUE)
map(database = "state",regions = regions$ot , col = "purple", fill=T, add=TRUE)
dev.off()

library(ape)
tree <- bt[[1]]
meta <- strsplit(tree$tip.label, "_")
tipabb <- sapply(meta, "[[", 1L)
tree$tip.label <- tipabb
tipabb <- factor(tipabb)
levels(tipabb) <- levs
color <- factor(tipabb, levels=c("se", "mw", "ot"), labels=c("blue", "orange", "purple"))

png("raxml-tree.png", width=1000, height=1000)
plot(tree, type='fan', show.tip.label=TRUE, tip.color=as.character(color))
dev.off()

most_recent <- max(as.numeric(sapply(meta, "[[", 3)))
root_tip_dists <- penaltree:::get_nodeheights(tree_info[[1]]$btr, dist_from_root = TRUE)$tipheights
nms <- names(root_tip_dists)
tip_year <- tree_info[[1]]$td + most_recent

png("temporal-signal.png", width=1000, height=600, pointsize=20)
plot(root_tip_dists, tip_year[nms], ylab="Year",
     xlab="Root to tip distance (sub. / site)", bty="l")
dev.off()

png("regpath-data.png", width=1000, height=750, pointsize=30)
plot(sp$fit, xvar="lambda")
dev.off()

library(c060)
png("stabpath-data.png", width=1000, height=1000, pointsize=30)
plot(sp)
dev.off()

png("stabpath-sim.png", width=1000, height=1000, pointsize=30)
plot(sp_sim)
dev.off()
