#!/usr/bin/Rscript

library(ape)

trRegion <- function(x) {
    ind <- match(x, state.abb)
    reg <- state.region[ind]
    regF <- factor(reg,
                   levels=c("Northeast", "South", "North Central", "West"),
                   labels=c('A', 'T', 'C', 'G'))
    regF <- as.character(regF)
    regF[is.na(regF)] <- 'N'
    regF[regF=='A'] <- 'N'
    regF[x=='Mexico'] <- 'A' ## let mexico have northeast label since no northeast sequences
    regF
}

tree <- read.nexus('mcc.tree')
nms <- tree$tip.label
abs <- sapply(strsplit(nms, '_'), '[[', 1)

regDNA <- trRegion(abs)
regDNA <- as.list(regDNA)
names(regDNA) <- nms

write.dna(regDNA, file='loc.aln', format='fasta')
write.tree(tree, file='mcc.nh')
save(regDNA, file='regDNA.RData')

flows <- read.csv("shipment-flows-origins-on-rows-dests-on-columns.csv", row.names=1)
                                        #TCGA
regRow <- trRegion(rownames(flows))
regCol <- trRegion(colnames(flows))

alphOrd <- c('T', 'C', 'G', 'A')
pairs <- expand.grid(to=alphOrd, from=alphOrd)
test <- pairs$from != pairs$to
pairs <- pairs[test,]
pairs <- pairs[, c('from', 'to')]

aggFlow <- function(from, to){
    rowTest <- regRow==from
    colTest <- regCol==to
    log10(sum(flows[rowTest, colTest]) + 1)
}
pairFlows <- mapply(aggFlow, from=pairs$from, to=pairs$to)
Z <- cbind("(Intercept)"=1, pairFlows)
cat(t(Z), file='designMat2')
