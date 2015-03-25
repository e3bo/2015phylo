library(ape)
library(boot)
library(ggplot2)
library(rphastRegression)

tmpf <- function(){
    tmNames <- system("grep ^TreeLikelihood beast/run1/beast-stdout | cut -d\'(\' -f2 | cut -d\')\' -f1 | cut -d\'-\' -f1", inter=TRUE)
    uniquePatterns <- system("grep \"unique pattern count\" beast/run1/beast-stdout | cut -d\' \' -f7", inter=TRUE)
    names(uniquePatterns) <- tmNames
    uniquePatterns
}
uniquePatterns <- tmpf()

tnames <- c('nonsIndel-aligned.fasta-gb', 'sIndel-aligned.fasta-gb')
tfiles <- paste0(names(uniquePatterns), '-aligned.fasta-gb.sample.trees')
trees <- lapply(tfiles, read.nexus)
# uncomment to speed up
# trees <- list(trees[[1]][1:3], trees[[2]][1:3])
flows <- read.csv("shipment-flows-origins-on-rows-dests-on-columns.csv", row.names=1)

nms <- lapply(trees, attr, which='TipLabel')

tmpf <- function(x) {
    x <- strsplit(x, '_')
    sapply(x, '[', 1)
}
abs <- lapply(nms, tmpf)

absF <- factor(unlist(abs))
levs <- levels(absF)
n <- length(unlist(levs))
alph <- LETTERS[seq_len(n)]

absFt <- lapply(abs, factor, levels=levs)
absFt <- lapply(absFt, 'levels<-', value=alph)
alph <- paste(alph, collapse='')
seqs <- lapply(absFt, as.character)
tmpf <- function(x, y) {
    msa(seqs=x, names=y, alphabet=alph)
}
pedvMSA <- mapply(tmpf, seqs, nms, SIMPLIFY=FALSE)

tmpf <- function(x) {
    unname(lapply(x, write.tree))
}
treeChar <- lapply(trees, tmpf)

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

tmpf <- function(x) {
    lapply(x, tm, subst.mod="UNREST", alphabet=alph, backgd=bg)
}
mods <- lapply(treeChar, tmpf)

assignElement <- function(x, nam, val) {
    tmpf <- function(x) {
        x[[nam]] <- val
        x
    }
    lapply(x, tmpf)
}
mods <- lapply(mods, assignElement, nam='rate.matrix', val=matrix(NA, nrow=n, ncol=n))

M <- list()
M[['sym']] <- lapply(mods, assignElement, nam='design.matrix', val=designMatrixSym)
M[['asym']] <- lapply(mods, assignElement, nam='design.matrix', val=designMatrixAsym)

getRateMatrix <- function(design.matrix, w){
    eta <- exp(design.matrix %*% w)
    pos <- 1
    rate.matrix <- matrix(nrow=n, ncol=n)
    for(i in seq_len(n)){
        rowSum <- 0
        for(j in seq_len(n)){            
            if (i != j) {
                rate.matrix[i,j] <- eta[pos]
                rowSum <- rowSum + eta[pos]
                pos <- pos + 1
            }            
        }
        rate.matrix[i,i] <- -rowSum        
    }
    rate.matrix
}

obj <- function(w, msal=pedvMSA, tmlol=M[['asym']]){
    design.matrix <- tmlol[[1]][[1]][['design.matrix']]
    rate.matrix <- getRateMatrix(design.matrix, w)
    tmlol <- lapply(tmlol, assignElement, nam='rate.matrix', val=rate.matrix)
    tmpf <- function(x, tmlist) {
        treeLogLikGivenMSA <- function(tm) likelihood.msa(x=x, tm=tm)
        sapply(tmlist, treeLogLikGivenMSA)
    }
    ll <- mapply(tmpf, msal, tmlol)
    mean(rowSums(ll))
}

ans <- list()
system.time(ans[['asym']] <- optim.rphast(obj, c(.001,.002), lower=c(-4,-2), upper=c(2,2)))
system.time(ans[['sym']] <- optim.rphast(obj, tmlol=M[['sym']], c(-1,.4), lower=c(-4,-2), upper=c(2,2)))
objNull <- function(x) obj(c(x, 0))
system.time(ans[['null']] <- optimize(objNull, interval=c(-4,2), maximum=TRUE))

(Dsym <- -2*ans[['null']]$objective + 2*-ans[['sym']]$value)
(Dasym <- -2*ans[['null']]$objective + 2*-ans[['asym']]$value)

#' A chi squared test supports rejection of the null for the directed model
pchisq(q=c(Dsym, Dasym), df=1, lower.tail=FALSE)

#' Parametric bootstrap for confidence intervals and bias

simPars <- ans[['asym']]$par

ran.gen <- function(data, pars, tmlol=M[['asym']], nsim=1){
    ntrees <- sapply(tmlol, length)
    ind <- sapply(ntrees, sample.int, size=1)
    tmpf <- function(x, y) x[[y]]
    tmsel <- mapply(tmpf, tmlol, ind, SIMPLIFY=FALSE)
    tmpf <- function(x) {
        design.matrix <- x[['design.matrix']]
        x[['rate.matrix']] <- getRateMatrix(design.matrix, pars)
        simulate.msa(object=x, nsim=nsim)
    }
    lapply(tmsel, tmpf)
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
est['point', ] <- ans[['asym']]$par
est['width', ] <- est['upper', ] - est['lower',]
est['plotMin', ] <- est['point', ] - est['width', ]/2*1.5
est['plotMax', ] <- est['point', ] + est['width', ]/2*1.5 

xlim <- est[c('plotMin', 'plotMax'), 'int']
ylim <- est[c('plotMin', 'plotMax'), 'flo']

D <- expand.grid(Intercept=seq(from=xlim[1], to=xlim[2], length.out=41),
                 Flows=seq(from=ylim[1], to=ylim[2], length.out=31))
D[, 'Expected\nlog likelihood'] <- apply(D, 1, obj)

estCol <- "#CC6677"
theme_set(theme_classic())
g <- ggplot(data=D, aes(x=Intercept, y=Flows, z=`Expected\nlog likelihood`))
g <- g + geom_tile(aes(fill=`Expected\nlog likelihood`))
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

nsim <- 1e3
simMsa <- ran.gen(data=NA, pars=simPars, nsim=nsim)

tmpf <- function(x) {
    x <- floor(10^x)
    cols <- seq_len(x)
    tmpff <- function(xx) {
        xx[, cols]
    }
    msa <- lapply(simMsa, tmpff)
    ans <- optim.rphast(obj, params=simPars, lower=c(-5,-5), upper=c(2,2), msa=msa)
    ans$par
}

inds <- seq(from=0, to=log10(nsim), length.out=20)
system.time(parSeq <- sapply(inds, tmpf))

par(mfrow=c(2,1))
par(mar=c(5,5,4,2) + .1)
plot(c(inds, inds[1]), c(parSeq[1, ], simPars[1]), type='n',
     xlab='log10(Fold increase in information)',
     ylab='Intercept\n( log {interstate movements} / year)')
points(inds, parSeq[1, ], type='b')
abline(h=simPars[1])
plot(c(inds, inds[1]), c(parSeq[2, ], simPars[2]), type='n',
     xlab='log10(Fold increase in information)',
     ylab='Flow effect\n(log{rate multiplier} / log {flow})')
points(inds, parSeq[2, ], type='b')
abline(h=simPars[2])

save.image('mle.RData')
