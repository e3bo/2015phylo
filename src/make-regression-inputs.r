#!/usr/bin/Rscript

library(ape)

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

n <- length(alph)
modFile <- 'init.mod'
cat("ALPHABET:", alph, "\n", file=modFile)
cat("ORDER: 0", "\n", file=modFile, append=TRUE)
cat("SUBST_MOD: UNREST", "\n", file=modFile, append=TRUE)
cat("TRAINING_LNL: 0", "\n", file=modFile, append=TRUE) ## Not sure all of these tags are needed
cat("BACKGROUND:", rep(1, times=n)/n, "\n", file=modFile, append=TRUE)
cat("RATE_MAT:", "\n", file=modFile, append=TRUE)
M <- matrix(1/(n-1), nrow=n, ncol=n)
diag(M) <- -1L
for(i in 1:n){
    cat(M[i,], sep="\t", file=modFile, append=TRUE)
    cat("\n", file=modFile, append=TRUE)
}
cat("TREE:", write.tree(tree), "\n", file=modFile, append=TRUE)


