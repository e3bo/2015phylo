#!/usr/bin/Rscript

library(ape)

trRegion <- function(x) {
    regF <- factor(x)
    alphabetSize <- length(levels(regF))
    ind <- seq_len(alphabetSize)
    levels(regF) <- LETTERS[ind]
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
flows <- read.csv("shipment-flows-origins-on-rows-dests-on-columns.csv", row.names=1)
write.tree(tree, file='mcc.nh')

nms <- tree$tip.label
abs <- sapply(strsplit(nms, '_'), '[[', 1)
absF <- factor(abs)
levs <- levels(absF)

alph <- LETTERS[seq_along(levs)]
levels(absF) <- alph

regDNA <- lapply(absF, as.character)
names(regDNA) <- nms

write.dna(regDNA, file='sim.fasta', format='fasta')
save(regDNA, file='regDNA.RData')

pairs <- expand.grid(to=levs, from=levs)
test <- pairs$from != pairs$to
pairs <- pairs[test,]
pairs <- pairs[, c('from', 'to')]

test <- abs %in% state.abb
pmf <- table(abs[test])
pmf <- pmf/sum(pmf)
pnms <- names(pmf)

usaRow <- colSums(flows[pnms, ] * as.numeric(pmf))
usaCol <- rowSums(t(t(flows[, pnms]) * as.numeric(pmf)))

M <- cbind(flows, 'USA'=usaCol)
M <- rbind(M, 'USA'=c(usaRow, 0))

aggFlow <- function(from, to, sym=TRUE){
    tot <- M[from, to]
    if(sym){
        tot <- tot + M[to, from]
    }
    log10(tot + 1)
}

pairFlows <- mapply(aggFlow, from=pairs$from, to=pairs$to)
Z <- cbind("(Intercept)"=1, pairFlows)
cat(t(Z), file='designMat2')




