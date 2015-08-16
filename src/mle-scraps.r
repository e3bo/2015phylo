nll <- function(x) -obj(x, tmlol=M[["big"]])
migsPerTime <- getInit()
par <- c(migsPerTime, rep(0, nc+1))
system.time(my.ans08 <- dtnet(F=nll, par=par, r=0.01, maxIter=100, a=0.1, tol=0.001, verbose=TRUE, debug=TRUE,
                             nlambda=100, log10LambdaRange=2, relStart=0.1, beta=0.5, mubar=1, alpha=0.8))

all(sapply(my.ans08, '[[', 'convergence')=='yes')
(lambda <- sapply(my.ans08, '[[', 'lambda'))
(kvec <- sapply(my.ans08, '[[', 'k'))
(sapply(my.ans08, '[[', 'nsg'))
(sapply(my.ans08, '[[', 'mu'))

my.path <- sapply(my.ans08, '[[', 'par')
matplot(lambda, t(my.path), type='l', log='x')
matplot(lambda, t(my.path[-1,]), type='l', log='x')
dev.off()

hc <- getClusters(migsPerTime=exp(migsPerTime))
clusts <- cutree(hc, k=10)

## Check that clusters are spread out on tree tips
t1 <- read.tree(text=M[['asym']][[1]][[1]]$tree)
t2 <- read.tree(text=M[['asym']][[2]][[1]]$tree)
trs <- c(t1, t2)
plot(trs[[1]], tip.color=clusts[trs[[1]]$tip.label])
plot(trs[[2]], tip.color=clusts[trs[[2]]$tip.label])

plot(read.tree(text=prune.tree(M[['asym']][[2]][[1]]$tree, names(clusts)[clusts==2])))

filterTips <- function(tmlol, keepers){
    prune <- function(xx){
        tr <- xx$tree
        tr <- prune.tree(tr, keepers, all.but=TRUE)
        xx$tree <- tr
        xx
    }
    tmpf <- function(x) {
        lapply(x, prune)
    }
    lapply(tmlol, tmpf)
}

nms <- names(clusts)[clusts!=2]
fold <- filterTips(M[['big']], nms)

filterSeqs <- function(msal, keepers){
    tmpf <- function(x){
        test <- x$names %in% keepers
        x[test, ]
    }
    lapply(msal, tmpf)
}
msaFold <- filterSeqs(pedvMSA, keepers=nms)
foldnll <- function(x) -obj(w=x, msal=msaFold, tmlol=fold)

system.time(fold.ans <- dtnet(F=foldnll, par=par, r=0.01, maxIter=100, a=0.1, tol=0.001, verbose=TRUE, debug=TRUE,
                             nlambda=100, log10LambdaRange=2, relStart=0.1, beta=0.99, mubar=1))

all(sapply(fold.ans, '[[', 'convergence')=='yes')
(foldlambda <- sapply(fold.ans, '[[', 'lambda'))
(foldkvec <- sapply(fold.ans, '[[', 'k'))
(sapply(fold.ans, '[[', 'nsg'))
(sapply(fold.ans, '[[', 'mu'))

fold.path <- sapply(fold.ans, '[[', 'par')
matplot(foldlambda, t(fold.path), type='l', log='x')
matplot(foldlambda, t(fold.path[-1,]), type='l', log='x')
dev.off()

foldnll <- sapply(fold.ans, '[[', 'F')
foldprednll <- apply(fold.path, 2, nll)
plot(foldlambda, -foldprednll + foldnll)
fold.path[,which.max(-foldprednll + foldnll)]

cv.dtnet <- function(msal, tmlol, nfolds){
    stopifnot(nfolds >2)
    totalNll <- function(x) -obj(x, tmlol=tmlol, msal=msal)
    migsPerTime <- getInit(msal=msal, tmlol=tmlol)
    nc <- ncol(tmlol[[1]][[1]]$design.matrix)
    parInit <- c(migsPerTime, rep(0, nc))
    hc <- getClusters(migsPerTime=exp(migsPerTime), tmlol=tmlol)
    clusts <- cutree(hc, k=nfolds)
    totAns <- dtnet(F=totalNll, par=parInit, r=0.01, maxIter=100,
                     a=0.1, tol=0.001, verbose=TRUE, debug=TRUE,
                     nlambda=100, log10LambdaRange=2, relStart=0.1,
                     beta=0.99, mubar=1, lower.limits=-5, upper.limits=5)
    lambda <- sapply(totAns, '[[', 'lambda')
    tmpf <- function(foldid){
        test <- clusts != foldid
        keepers <- names(clusts)[test]
        msalF <- filterSeqs(msal, keepers=keepers)
        tmlolF <- filterTips(tmlol, nms)
        foldNll <- function(x) -obj(w=x, msal=msalF, tmlol=tmlolF)
        foldAns <- dtnet(F=foldNll, par=parInit, r=0.01, maxIter=100,
                           a=0.1, tol=0.001, verbose=TRUE, debug=TRUE,
                          lambda=lambda, beta=0.99, mubar=1,
                          lower.limits=-5, upper.limits=5)
        fittedNll <- sapply(foldAns, '[[', 'F')
        conv <- sapply(foldAns, '[[', 'convergence')
        stopifnot(all(conv=='yes'))
        k <- sapply(foldAns, '[[', 'k')
        foldPath <- sapply(foldAns, '[[', 'par')
        predNll <- apply(foldPath, 2, totalNll)
        cv <- -predNll - (-fittedNll)
        data.frame(foldid=foldid, lambda=lambda, fittedNll=fittedNll, predNll=predNll, cv=cv, k=k)
    }
    cv <- lapply(seq_len(nfolds), tmpf)
    list(totAns, cv)
}

cvasym <- cv.dtnet(msal=pedvMSA, tmlol=M[['asym']], nfolds=10)
cvs <- sapply(cvasym[[2]], '[[', 'cv')
EstPredLL <- rowMeans(cvs)
(PredLLSE <- apply(cvs, 1, sd)/sqrt(ncol(cvs)))
plot(cvasym[[2]][[1]][, 'lambda'], EstPredLL, log='x')

totpath <- sapply(cvasym[[1]], '[[', 'par')

#cvasym <- do.call(rbind, cvasym)



#' ## Check for simulation bias
#'

ran.gen <- function(data, pars, tmlol=M[['asym']], nsim=1){
    ntrees <- sapply(trees, length)
    ind <- sapply(ntrees, sample.int, size=1)
    tmpf <- function(x, tr, y) {
        ret <- x[[1]]
        ret$tree <- write.tree(tr[y])
        ret
    }
    tmsel <- mapply(tmpf, tmlol, trees, ind, SIMPLIFY=FALSE)
    tmpf <- function(x) {
        design.matrix <- x[['design.matrix']]
        x[['rate.matrix']] <- getRateMatrix(design.matrix, pars)
        simulate.msa(object=x, nsim=nsim)
    }
    lapply(tmsel, tmpf)
}

get.score.stat <- function(data, pars) {
    grad(obj, x=pars, method='simple', msal=data)
}

simulation.bias.diagnostic <- function(pars, R=1e3) {
    bs.score <- boot(data=pedvMSA, get.score.stat, R=R, sim='parametric',
                     ran.gen=ran.gen, mle=pars, parallel='multicore',
                     ncpus=parallel::detectCores(), pars=pars)
    print(bs.score)
    hats <- colMeans(bs.score$t)
    V <- var(bs.score$t)
    Q <- hats %*% solve(V) %*% hats
    print(Q)
    1 - pchisq(q=Q, df=length(pars))
}

simPars <- lapply(ans, '[[', 'par')
simPars$null <- c(ans$null$maximum, 0)

(pvals <- sapply(simPars, simulation.bias.diagnostic, R=1e3))
stopifnot(pvals > 0.1)

#' ## Likelihood ratio test

(Dsym <- -2*ans[['null']]$objective + 2*-ans[['sym']]$value)
(Dasym <- -2*ans[['null']]$objective + 2*-ans[['asym']]$value)

#' A chi squared test supports rejection of the null for the directed model
pchisq(q=c(Dsym, Dasym), df=1, lower.tail=FALSE)

#' ## Parametric bootstrap for confidence intervals and bias

get.param.stat <- function(data, pars) {
    ans <- optim.rphast(obj, params=pars, lower=c(-5,-5), upper=c(2,2), msal=data)
    ans$par
}

bsParamR <- 1e4
system.time(bs <- boot(data=pedvMSA, get.param.stat, R=bsParamR, sim='parametric',
                       ran.gen=ran.gen, mle=simPars$asym, pars=simPars$asym, parallel='multicore',
                       ncpus=parallel::detectCores()))

# Bootstrap bias and standard error estimates

bs

# The distribution of estimates appears close to normal, without any discontinuities
plot(bs, index=1)
plot(bs, index=2)

(ciInt <- boot.ci(bs, index=1, type=c('perc', 'norm')))
(ciFlo <- boot.ci(bs, index=2, type=c('perc', 'norm')))

#' ##  Plot of the point estimates, confidence intervals, and likelihood surface

est <- data.frame(int=ciInt$percent[4:5], flo=ciFlo$percent[4:5], row.names=c('lower', 'upper'))
est['point', ] <- ans[['asym']]$par
est['width', ] <- est['upper', ] - est['lower',]
est['plotMin', ] <- est['point', ] - est['width', ]/2*1.5
est['plotMax', ] <- est['point', ] + est['width', ]/2*1.5

xlim <- est[c('plotMin', 'plotMax'), 'int']
ylim <- est[c('plotMin', 'plotMax'), 'flo']

D <- expand.grid(Intercept=seq(from=xlim[1], to=xlim[2], length.out=41),
                 Flows=seq(from=ylim[1], to=ylim[2], length.out=31))
D[, 'Simulated\nlog likelihood'] <- apply(D, 1, obj)

estCol <- "#CC6677"
theme_set(theme_classic())
g <- ggplot(data=D, aes(x=Intercept, y=Flows, z=`Simulated\nlog likelihood`))
g <- g + geom_tile(aes(fill=`Simulated\nlog likelihood`))
g <- g + stat_contour(colour="#DDCC77", alpha=0.5)
g <- g + geom_point(x=est['point', 'int'], y=est['point', 'flo'], size=5, colour=estCol)
g <- g + xlab('\nScale') + ylab('Flow effect\n')
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

##' ## Check consistency

nsim <- 1e3
## This simulation method is not correct (has been outdated by other developments), and so simulation bias will be evident
simMsa <- ran.gen(data=NA, pars=simPars$asym, nsim=nsim)

tmpf <- function(x) {
    x <- floor(10^x)
    cols <- seq_len(x)
    tmpff <- function(xx) {
        xx[, cols]
    }
    msa <- lapply(simMsa, tmpff)
    ans <- optim.rphast(obj, params=simPars$asym, lower=c(-5,-5), upper=c(2,2), msa=msa)
    ans$par
}

inds <- seq(from=0, to=log10(nsim), length.out=20)
system.time(parSeq <- sapply(inds, tmpf))

par(mfrow=c(2,1))
par(mar=c(5,5,4,2) + .1)
plot(c(inds, inds[1]), c(parSeq[1, ], simPars$asym[1]), type='n',
     xlab='log10(Fold increase in information)',
     ylab='Intercept\n( log {interstate movements} / year)')
points(inds, parSeq[1, ], type='b')
abline(h=simPars$asym[1])
plot(c(inds, inds[1]), c(parSeq[2, ], simPars$asym[2]), type='n',
     xlab='log10(Fold increase in information)',
     ylab='Flow effect\n(log{rate multiplier} / log {flow})')
points(inds, parSeq[2, ], type='b')
abline(h=simPars$asym[2])
