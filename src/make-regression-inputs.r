#!/usr/bin/Rscript

library(ape)

tree <- read.nexus('mcc.tree')
nms <- tree$tip.label
abs <- sapply(strsplit(nms, '_'), '[[', 1)
ind <- match(abs, state.abb)
reg <- state.region[ind]

regDNA <- factor(reg,
                 levels=c("Northeast", "South", "North Central", "West"),
                 labels=c('A', 'T', 'C', 'G'))
regDNA <- as.character(regDNA)
regDNA[is.na(regDNA)] <- 'N'
regDNA[abs=='Mexico'] <- 'A'
## let mexico have northeast label since no northeast sequences
regDNA <- as.list(regDNA)
names(regDNA) <- nms
write.dna(regDNA, file='loc.aln', format='fasta')
write.tree(tree, file='mcc.nh')
