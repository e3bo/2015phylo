#!/usr/bin/Rscript

library(rphast)

tree <- read.tree('mcc.nh')
load('regDNA.RData')

locMsa <- msa(regDNA, names(regDNA), alphabet='ACTGN')
treeChar <- write.tree(tree)

mod <- phyloFit(locMsa, tree=treeChar, subst.mod='UNREST', no.opt='branches', ninf.sites=1)
mod$tree <- treeChar

mod2 <- phyloFit(locMsa, init.mod=mod, subst.mod='UNREST', no.opt='branches', ninf.sites=1)

pdf('unrest-rate-estimates.pdf')
plot(mod2)
dev.off()
