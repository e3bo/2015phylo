library(ape)
library(rphastRegression)
library(ggplot2)

tree <- read.tree('mcc.nh')
flows <- read.csv("shipment-flows-origins-on-rows-dests-on-columns.csv", row.names=1)

nms <- tree$tip.label
abs <- sapply(strsplit(nms, '_'), '[[', 1)
absF <- factor(abs)
levs <- levels(absF)

n <- length(levs)
alph <- LETTERS[seq_len(n)]
levels(absF) <- alph
alph <- paste(alph, collapse='')
seqs <- lapply(absF, as.character)
pedvMSA <- msa(seqs=seqs, names=nms, alphabet=alph)
treeChar <- write.tree(tree)

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
designMatrix <- cbind("(Intercept)"=1, pairFlows)

bg <- pmf[levs]
bg[is.na(bg)] <- 0

mod <- tm(treeChar, "UNREST", alphabet=alph, backgd=bg)
mod$rate.matrix <- matrix(NA, nrow=n, ncol=n)
mod$design.matrix <- designMatrix

obj <- function(w){
    eta <- exp(mod$design.matrix %*% w)
    pos <- 1
    for(i in seq_len(n)){
        rowSum <- 0
        for(j in seq_len(n)){            
            if (i != j) {
                mod$rate.matrix[i,j] <- eta[pos]
                rowSum <- rowSum + eta[pos]
                pos <- pos + 1
            }            
        }
        mod$rate.matrix[i,i] <- -rowSum        
    }
    mod$rate.matrix
    likelihood.msa(pedvMSA, tm=mod)
}
obj(c(.001,.002))

ans <- optim.rphast(obj, c(.001,.002), lower=c(-2,-2), upper=c(2,2))
D <- expand.grid(x=seq(from=-1.5, to=-.5, length.out=20),
                 y=seq(from=-1, to=0.1, length.out=20))
D$z <- apply(D, 1, obj)

g <- ggplot(data=D, aes(x=x, y=y, z=z))
g <- g + geom_tile(aes(fill=z))
g <- g + stat_contour() + geom_point(x=ans$par[1], y=ans$par[2])
ggsave(g, 'mle.pdf')

q('no')

root.seq <- NucleotideSequence(length=10)
p<-GTR(rate.params=list("a"=1, "b"=1, "c"=3,"d"=2, "e"=1, "f"=1), base.freqs=c(1,1,1,1)/4)
attachProcess(root.seq,p)
sampleStates(root.seq)

sim <- PhyloSim()
sim$phylo <- tree
sim$rootSeq <- root.seq
Simulate(sim)
saveAlignment(sim,file="sim10.fas", skip.internal=TRUE)

