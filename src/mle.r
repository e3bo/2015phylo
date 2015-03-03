library(ape)
library(rphastRegression)
library(ggplot2)

tree <- read.nexus('mcc.tree')
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

aggFlow <- function(from, to, sym=TRUE, M=flows){
    tot <- M[from, to]
    if(sym){
        tot <- tot + M[to, from]
    }
    log10(tot + 1)
}

pairFlows <- mapply(aggFlow, from=as.character(pairs$from), to=as.character(pairs$to))
designMatrix <- cbind("(Intercept)"=1, pairFlows)

bg <- rep(1, n)/n

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
    likelihood.msa(pedvMSA, tm=mod)
}

ans <- optim.rphast(obj, c(.001,.002), lower=c(-4,-2), upper=c(2,2))
objNull <- function(x) obj(c(x, 0))
ansNull <- optimize(objNull, interval=c(-4,2), maximum=TRUE)
D <- -2*ansNull$objective + 2*-ans$value

#' A chi squared test does not allow for rejection of the null hypothesis
pchisq(q=D, df=1, lower.tail=FALSE)

D <- expand.grid(Intercept=seq(from=-5.5, to=1.5, length.out=41),
                 Flows=seq(from=-1., to=1., length.out=31))
D$logLikelihood <- apply(D, 1, obj)

#' Though the point estimated is positive, there is a strong correlation with the intercept
theme_set(theme_classic())
g <- ggplot(data=D, aes(x=Intercept, y=Flows, z=logLikelihood))
g <- g + geom_tile(aes(fill=logLikelihood))
g <- g + stat_contour() + geom_point(x=ans$par[1], y=ans$par[2])
g <- g + xlab('Intercept') + ylab('Flow effect')
g

#' The log likelihood is not convex far from the optimum
plot(logLikelihood~Flows, data=D[abs(D$Intercept-0.1) < .001,], type='l')

if(FALSE) {

root.seq <- NucleotideSequence(length=10)
p<-GTR(rate.params=list("a"=1, "b"=1, "c"=3,"d"=2, "e"=1, "f"=1), base.freqs=c(1,1,1,1)/4)
attachProcess(root.seq,p)
sampleStates(root.seq)

sim <- PhyloSim()
sim$phylo <- tree
sim$rootSeq <- root.seq
Simulate(sim)
saveAlignment(sim,file="sim10.fas", skip.internal=TRUE)
}
