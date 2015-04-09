library(ape)
library(boot)
library(ggplot2)
library(numDeriv)
library(rphastRegression)

tmpf <- function(){
    tmNames <- system("grep ^TreeLikelihood beast/run1/beast-stdout | cut -d\'(\' -f2 | cut -d\')\' -f1 | cut -d\'-\' -f1", inter=TRUE)
    uniquePatterns <- system("grep \"unique pattern count\" beast/run1/beast-stdout | cut -d\' \' -f7", inter=TRUE)
    names(uniquePatterns) <- tmNames
    uniquePatterns
}
uniquePatterns <- tmpf()

tnames <- c('nonsIndel', 'sIndel')
tfiles <- paste0(tnames, '-aligned.fasta-gb.combined.trees')
trees <- lapply(tfiles, read.nexus)

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

simulationR <- 2
tmpf <- function(x) {
    n <- length(x)
    ind <- sample.int(n, size=simulationR, replace=TRUE)
    x <- x[ind]
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
designMatrixSym <- data.matrix(pairFlows)

pairFlows <- mapply(aggFlow, from=as.character(pairs$from), to=as.character(pairs$to), sym=FALSE)
designMatrixAsym <- data.matrix(pairFlows)

nr <- nrow(designMatrixAsym)
nc <- 10
mNoise <-matrix(runif(nr*nc), nrow=nr)
designMatrixBigger <- cbind(designMatrixAsym, mNoise)
designMatrixBigger <- scale(designMatrixBigger)

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
M[['big']] <- lapply(mods, assignElement, nam='design.matrix', val=designMatrixBigger)

##' ## Fit models

getRateMatrix <- function(design.matrix, w){
    scale <- exp(w[1])
    effects <- w[-1]
    eta <- exp(design.matrix %*% effects)
    eta <- eta / mean(eta) * scale
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
    ll <- rowSums(ll)
    scale <- max(ll)
    ll <- ll - scale
    probs <- exp(ll)
    log(mean(probs)) + scale
}

getInit <- function(msal=pedvMSA, tmlol=M[['asym']]){
    tmpf <- function(x, y){
        dloc <- outer(y$seq, y$seq, '==')
        colnames(dloc) <- names(y)
        rownames(dloc) <- names(y)
        tmpff <- function(xx) {
            phy <- read.tree(text=xx$tree)
            dphy <- cophenetic(phy)
            dphy <- dphy[rownames(dloc), colnames(dloc)]
            rdist <- dloc / dphy
            mean(rdist[upper.tri(rdist)])
        }
        lapply(x, tmpff)
    }
    res <- mapply(tmpf, x=tmlol, y=msal, SIMPLIFY=FALSE)
    tmpf <- function(x) mean(unlist(x))
    res <- sapply(res, tmpf)
    tmpf <- function(x) length(x$seq)
    nseq <- sapply(pedvMSA, tmpf)
    w <- nseq*(nseq - 1)
    res <- weighted.mean(res, w=w)
    log(res)
}

getClusters <- function(tmlol=M[['asym']], migsPerTime=1){
    tmpf <- function(x, y){
        tmpff <- function(xx) {
            phy <- read.tree(text=xx$tree)
            cophenetic(phy)
        }
        res <- lapply(x, tmpff)
        res <- Reduce('+', res)/length(res)
        -res
    }
    res <- lapply(tmlol, tmpf)
    joinMats <- function(A, B){
        od <- matrix(-Inf, nrow=nrow(B), ncol=ncol(A))
        rownames(od) <- rownames(B)
        colnames(od) <- colnames(A)
        left <- rbind(A, od)
        od <- matrix(-Inf, nrow=nrow(A), ncol=ncol(B))
        rownames(od) <- rownames(A)
        colnames(od) <- colnames(B)
        right <- rbind(od, B)
        cbind(left, right)
    }
    res <- Reduce(joinMats, res)
    res <- exp(res*migsPerTime)
    d <- as.dist(res)
    hc <- hclust(d, method='complete')
    hc
}

my.opt <- function(F, par, maxIter=100, tol=1e-2, a=0.1, r=0.01, upper=2,
                   lower=-2, relStart=1, nlambda=1, log10LambdaRange=2,
                   mubar=1, beta=.9, verbose=FALSE, debug=TRUE, initFactor=10){
    niter <- 0
    dim <- length(par)
    parInds <- 1:dim
    I <- diag(nrow=dim)
    gF <- grad(F, x=par, method='simple')
    mu <- mubar
    stopifnot(beta >0, beta <1)
    G <- diag(initFactor * abs(gF), ncol=dim)
    lstart <- max(abs(gF))
    loglstart <- log10(lstart) + relStart
    loglend <- loglstart - log10LambdaRange
    logstep <- (loglend - loglstart)/nlambda
    lscaler <- 10^logstep
    lambda <- c(0, rep(10^loglstart, dim-1))
    lower <- rep(lower, dim)
    upper <- rep(upper,dim)
    res <- list()
    fsg <- function(p, g, l, h) {
        if(p==0) {
            max(abs(g) - l, 0)
        } else if (p > 0){
            (-g - l)/h
        } else {
            (-g + l)/h
        }
    }
    for (i in seq_len(nlambda)){
        k <- 0
        F1 <- F(par) + abs(par) %*% lambda
        H <- I/(2*mu) + G
        sg <- mapply(fsg, p=par, g=gF, l=lambda, h=diag(H))
        nsg <- max(abs(sg))
        while (nsg > tol && k < maxIter){
            d <- numeric(dim)
            dlist <- list()
            fmlist <- list()
            inactive <- par != 0 | sg != 0 ## inactive means not actively fixed to zero in solution
            nInactive <- sum(inactive)
            ndesc <- (1 + floor(k * a)) * nInactive
            for (nd in 1:ndesc){
                j <- sample.int(n=nInactive, size=1)
                j <- parInds[inactive][j]
                Hd <- H %*% d
                gr <- gF[j] + Hd[j]
                if (par[j] + d[j] > 0 || (par[j] + d[j] == 0 & -gr > 0)){
                    z <- (-gr - lambda[j])/H[j,j]
                    if (par[j] + d[j] + z < 0){
                        d[j] <- -par[j]
                    } else {
                        d[j] <- d[j] + z
                    }
                } else {
                    z <- (-gr + lambda[j])/H[j,j]
                    if (par[j] + d[j] + z > 0){
                        d[j] <- -par[j]
                    } else {
                        d[j] <- d[j] + z
                    }
                }
                dlist[[nd]] <- d
                fmlist[[nd]] <- d %*% (H/2) %*% d + gF %*% d + F1 - abs(par) %*% lambda + abs(par + d) %*% lambda
            }
            if(any(diff(c(F1, unlist(fmlist)))>1e-10)) browser()
            par2 <- par + d
            if (any(par2 > upper) || any(par2 < lower)){
                if(verbose) cat('backtracking: out of bounds', '\n')
                mu <- mu * beta
            } else {
                Fmod <- d %*% (H/2) %*% d + gF %*% d + F1 - abs(par) %*% lambda + abs(par2) %*% lambda
                if(Fmod  > F1 && !isTRUE(all.equal(Fmod, F1))) {
                    if (debug) browser()
                    if(verbose) cat('backtracking: poorly solved model', '\n')
                    mu <- mu * beta
                } else {
                    try(F2 <- F(par2) + abs(par2) %*% lambda)
                    if(inherits(F2, 'try-error')) browser()
                    if (F2 - F1 > r * (Fmod - F1)){
                        if(verbose) cat('backtracking: insufficient decrease', '\n')
                        mu <- mu * beta
                    } else {
                        gF2 <- grad(F, x=par2, method='simple')
                        y <- gF2 - gF
                        s <- d
                        ys <- y%*%s
                        if (ys > 0){
                            ## Hessian approximation update via BFGS via 8.19 in Nocedal and Wright
                            sColVec <- matrix(s, ncol=1)
                            yColVec <- matrix(y, ncol=1)
                            M <- G %*% sColVec %*% t(sColVec) %*% G
                            M <- M / as.numeric(t(sColVec) %*% G %*% sColVec)
                            M <- G - M
                            G <- M + yColVec %*% t(yColVec) / as.numeric(t(yColVec) %*% sColVec)
                            if (any(diag(G) <= 0)) browser()
                        } else if(verbose){
                            cat('skipping Hessian update: ys <= 0', '\n')
                        }
                        ## update vars
                        k <- k + 1
                        par <- par2
                        gF <- gF2
                        dF <- F2 - F1
                        F1 <- F2
                        H <- I/(2*mu) + G
                        sg <- mapply(fsg, p=par, g=gF, l=lambda, h=diag(H))
                        nsg <- max(abs(sg))
                        if (debug && mu < 1e-8) browser()
                        if (verbose) {
                            cat('lambda[2]: ', lambda[2], '\n')
                            cat('k: ', k, '\n')
                            cat('F: ', as.numeric(F1), '\n')
                            cat('par: ', signif(par, 3), '\n')
                            cat('grad: ', signif(gF, 3), '\n')
                            cat('nsg: ', signif(nsg, 3), '\n')
                            cat('\n')
                        }
                    }
                }
            }
        }
        convergence <- ifelse(k == maxIter, 'no', 'yes')
        res[[i]] <- list(par=par, F=F2, k=k, gF=gF2, H=H, lambda=lambda, convergence=convergence, mu=mu, nsg=nsg, sg=sg)
        mu <- mubar
        lambda <- lambda * lscaler
    }
    res
}

nll <- function(x) -obj(x, tmlol=M[["big"]])
migsPerTime <- getInit()
par <- c(migsPerTime, rep(0, nc +1))
system.time(my.ans <- my.opt(F=nll, par=par, r=0.01, maxIter=100, a=0.1, tol=0.001, verbose=TRUE, debug=TRUE,
                             nlambda=100, log10LambdaRange=2, relStart=0.1, beta=0.99, mubar=1))

all(sapply(my.ans, '[[', 'convergence')=='yes')
(lambda <- sapply(my.ans, function(x) x$lambda[2]))
(kvec <- sapply(my.ans, '[[', 'k'))
(sapply(my.ans, '[[', 'nsg'))
(sapply(my.ans, '[[', 'mu'))

my.path <- sapply(my.ans, '[[', 'par')
matplot(lambda, t(my.path), type='l', log='x')
matplot(lambda, t(my.path[-1,]), type='l', log='x')
#dev.off()

hc <- getClusters(migsPerTime=exp(migsPerTime))
clusts <- cutree(hc, k=10)

## Check that clusters are spread out on tree tips
t1 <- read.tree(text=M[['asym']][[1]][[1]]$tree)
t2 <- read.tree(text=M[['asym']][[2]][[1]]$tree)
trs <- c(t1, t2)
plot(trs[[1]], tip.color=clusts[trs[[1]]$tip.label])
plot(trs[[2]], tip.color=clusts[trs[[2]]$tip.label])

plot(read.tree(text=prune.tree(M[['asym']][[2]][[1]]$tree, names(clusts)[clusts==2])))


ans <- list()
system.time(ans[['asym']] <- optim.rphast(obj, c(.001,.002), lower=c(-4,-2), upper=c(2,2)))
system.time(ans[['sym']] <- optim.rphast(obj, tmlol=M[['sym']], c(-1,.4), lower=c(-4,-2), upper=c(2,2)))
objNull <- function(x) obj(c(x, 0))
system.time(ans[['null']] <- optimize(objNull, interval=c(-4,2), maximum=TRUE))

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

save.image('mle.RData')
