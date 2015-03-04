library(ape)
library(rphastRegression)
library(ggplot2)
library(phylosim)

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

pairFlows <- mapply(aggFlow, from=as.character(pairs$from), to=as.character(pairs$to), sym=TRUE)
designMatrixSym <- cbind("(Intercept)"=1, pairFlows)

pairFlows <- mapply(aggFlow, from=as.character(pairs$from), to=as.character(pairs$to), sym=FALSE)
designMatrixAsym <- cbind("(Intercept)"=1, pairFlows)

bg <- rep(1, n)/n

mod <- tm(treeChar, "UNREST", alphabet=alph, backgd=bg)
mod$rate.matrix <- matrix(NA, nrow=n, ncol=n)
modSym <- modAsym <- mod
modSym$design.matrix <- designMatrixSym
modAsym$design.matrix <- designMatrixAsym

setRateMatrix <- function(mod, w){
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
    mod
}

obj <- function(w, msa=pedvMSA, mod=modAsym){
    mod <- setRateMatrix(mod, w)
    likelihood.msa(msa, tm=mod)
}

ansAsym <- optim.rphast(obj, c(.001,.002), lower=c(-4,-2), upper=c(2,2))
ansSym <- optim.rphast(obj, mod=modSym, c(-1,.4), lower=c(-4,-2), upper=c(2,2))
objNull <- function(x) obj(c(x, 0))
ansNull <- optimize(objNull, interval=c(-4,2), maximum=TRUE)
Dsym <- -2*ansNull$objective + 2*-ansSym$value
Dasym <- -2*ansNull$objective + 2*-ansAsym$value

#' A chi squared test supports rejection of the null for the directed model
pchisq(q=c(Dsym, Dasym), df=1, lower.tail=FALSE)

D <- expand.grid(Intercept=seq(from=-5.5, to=1.5, length.out=41),
                 Flows=seq(from=-1., to=1., length.out=31))
D$logLikelihood <- apply(D, 1, obj)

#' There is a strong correlation with the intercept
theme_set(theme_classic())
g <- ggplot(data=D, aes(x=Intercept, y=Flows, z=logLikelihood))
g <- g + geom_tile(aes(fill=logLikelihood))
g <- g + stat_contour() + geom_point(x=ansAsym$par[1], y=ansAsym$par[2])
g <- g + xlab('Intercept') + ylab('Flow effect')
g

#' The log likelihood is not convex far from the optimum
plot(logLikelihood~Flows, data=D[abs(D$Intercept-0.1) < .001,], type='l')

#' Now check consistency with simulation
#'
alphav <- strsplit(alph, split=NULL)[[1]]
a <- Alphabet(alphav)

PSIM_FAST <- TRUE

simPars <- ansAsym$par
rm <- setRateMatrix(modAsym, simPars)$rate.matrix
expectedSubsPerTime <- -min(as.numeric(eigen(rm)$values))

rateNames <- outer(alphav,alphav, paste, sep="->")
rates <- as.numeric(rm)
names(rates) <- as.character(rateNames)
isDiag <- grepl("([A-Z])->\\1", names(rates))
rates <- rates[!isDiag]
gs <- GeneralSubstitution(name="geoSubs", alphabet=a, rate.list=as.list(rates))

root.seq <- Sequence(length=10, alphabets=list(a))
attachProcess(root.seq,gs)
sampleStates(root.seq)

sim <- PhyloSim()
sim$phylo <- tree
scaleTree(sim, expectedSubsPerTime)

sim$rootSeq <- root.seq
Simulate(sim)
saveAlignment(sim,file="sim.fasta", skip.internal=TRUE)

simMsa <- read.msa('sim.fasta', alphabet=alph)
ansSim <- optim.rphast(obj, params=simPars, lower=c(-4,-20), upper=c(2,2), msa=simMsa)

D$simLogLikelihood <- apply(D[, c('Intercept', 'Flows')], 1, obj, msa=simMsa)

#' The simulated data has a likelihood function that looks similar to that of the real data
theme_set(theme_classic())
g <- ggplot(data=D, aes(x=Intercept, y=Flows, z=simLogLikelihood))
g <- g + geom_tile(aes(fill=simLogLikelihood))
g <- g + stat_contour() + geom_point(x=ans$par[1], y=ans$par[2])
g <- g + xlab('Intercept') + ylab('Flow effect')
g

