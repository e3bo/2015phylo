library(ape)
library(boot)
library(ggplot2)
library(rphastRegression)


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
(Dsym <- -2*ansNull$objective + 2*-ansSym$value)
(Dasym <- -2*ansNull$objective + 2*-ansAsym$value)

#' A chi squared test supports rejection of the null for the directed model
pchisq(q=c(Dsym, Dasym), df=1, lower.tail=FALSE)


#' The log likelihood is not convex far from the optimum, which
#' illustrates the importance of starting the search in the convex
#' region around the optimum. I guess this is probably generally true
#' for optimization of likelihoods.

plot(logLikelihood~Flows, data=D[abs(D$Intercept-0.1) < .001,], type='l')

#' Parametric bootstrap for confidence intervals and bias

simPars <- ansAsym$par

ran.gen <- function(data, pars){
    msim <- setRateMatrix(modAsym, w=pars)
    simulate.msa(object=msim, nsim=1)
}

get.stat <- function(data) {
    ans <- optim.rphast(obj, params=simPars, lower=c(-5,-5), upper=c(2,2), msa=data)
    ans$par
}

system.time(bs <- boot(data=pedvMSA, get.stat, R=1e3, sim='parametric',
                       ran.gen=ran.gen, mle=simPars, parallel='multicore',
                       ncpus=parallel::detectCores()))

# The distribution of estimates appears close to normal, without any discontinuities
plot(bs, index=1)
plot(bs, index=2)

(ciInt <- boot.ci(bs, index=1, type=c('perc', 'norm')))
(ciFlo <- boot.ci(bs, index=2, type=c('perc', 'norm')))

#' Make a nice plot of the point estimates, confidence intervals, and likelihood surface

est <- data.frame(int=ciInt$percent[4:5], flo=ciFlo$percent[4:5], row.names=c('lower', 'upper'))
est['point', ] <- ansAsym$par
est['width', ] <- est['upper', ] - est['lower',]
est['plotMin', ] <- est['point', ] - est['width', ]/2*1.5
est['plotMax', ] <- est['point', ] + est['width', ]/2*1.5 

xlim <- est[c('plotMin', 'plotMax'), 'int']
ylim <- est[c('plotMin', 'plotMax'), 'flo']

D <- expand.grid(Intercept=seq(from=xlim[1], to=xlim[2], length.out=41),
                 Flows=seq(from=ylim[1], to=ylim[2], length.out=31))
D[, 'Log likelihood'] <- apply(D, 1, obj)

estCol <- "#CC6677"
theme_set(theme_classic())
g <- ggplot(data=D, aes(x=Intercept, y=Flows, z=`Log likelihood`))
g <- g + geom_tile(aes(fill=`Log likelihood`))
g <- g + stat_contour(colour="#DDCC77", alpha=0.5)
g <- g + geom_point(x=est['point', 'int'], y=est['point', 'flo'], size=5, colour=estCol)
g <- g + xlab('\nIntercept') + ylab('Flow effect\n')
g <- g + geom_segment(x=est['point', 'int'], xend=est['point', 'int'],
                      y=est['lower', 'flo'], yend=est['upper', 'flo'],
                      color=estCol, arrow=grid::arrow(angle=90, ends='both'),
                      lineend='round', size=1.25)
g <- g + geom_segment(x=est['lower', 'int'], xend=est['upper', 'int'],
                      y=est['point', 'flo'], yend=est['point', 'flo'],
                      colour=estCol, arrow=grid::arrow(angle=90, ends='both'),
                      lineend='round', size=1.25)
#g <- g + theme(legend.position='top')
ggsave('ll-surface.pdf', width=7, height=5, pointsize=12)

# Illustrate consistency

msim <- setRateMatrix(modAsym, w=simPars)
nsim <- 1e3
simMsa <- simulate.msa(object=msim, nsim=nsim)

tmpf <- function(x) {
    x <- floor(10^x)
    cols <- seq_len(x)
    msa <- simMsa[, cols]
    ans <- optim.rphast(obj, params=simPars, lower=c(-5,-5), upper=c(2,2), msa=msa)
    ans$par
}

inds <- seq(from=0, to=log10(nsim), length.out=20)
system.time(parSeq <- sapply(inds, tmpf))

par(mfrow=c(2,1))
par(mar=c(5,5,4,2) + .1)
plot(inds, parSeq[1, ], xlab='log10(Fold increase in information)',
     ylab='Intercept\n( log {interstate movements} / year)', type='b')
abline(h=simPars[1])
plot(inds, parSeq[2, ], xlab='log10(Fold increase in information)',
     ylab='Flow effect\n(log{rate multiplier} / log {flow})', type='b')
abline(h=simPars[2])

save.image('mle.RData')
