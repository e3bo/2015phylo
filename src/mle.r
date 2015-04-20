library(ape)
library(boot)
library(expm)
library(exptest)
library(ggplot2)
library(lubridate)
library(numDeriv)
library(rphastRegression)

#' ## Data loading

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

#' ## Tree-model creation

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

simulationR <- 1
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

##' ## Fit 1-parameter model

getRateMatrix <- function(design.matrix, w){
    scale <- exp(w[1])
    effects <- w[-1]
    stopifnot(ncol(design.matrix) == length(effects))
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
    if(!is.null(dim(ll))){
        ll <- rowSums(ll)
        scale <- max(ll)
        ll <- ll - scale
        probs <- exp(ll)
        ll <- log(mean(probs)) + scale
    } else if(length(ll) > 1){
        ll <- sum(ll)
    }
    ll
}

obj2 <- function(w, x, y){
    msal <- y$msal
    tmlol <- y$tmlol
    rate.matrix <- getRateMatrix(x, w)
    tmlol <- lapply(tmlol, assignElement, nam='rate.matrix', val=rate.matrix)
    tmpf <- function(xx, tmlist) {
        treeLogLikGivenMSA <- function(tm) likelihood.msa(x=xx, tm=tm)
        sapply(tmlist, treeLogLikGivenMSA)
    }
    ll <- mapply(tmpf, msal, tmlol)
    if(!is.null(dim(ll))){
        ll <- rowSums(ll)
        scale <- max(ll)
        ll <- ll - scale
        probs <- exp(ll)
        ll <- log(mean(probs)) + scale
    } else if(length(ll) > 1){
        ll <- sum(ll)
    }
    ll
}

objNull <- function(x) obj(c(x, 0))
system.time(oneParAns <- optimize(objNull, interval=c(-4,2), maximum=TRUE))

#' ## Fit 2-parameter models

ans <- list()
system.time(ans[['asym']] <- optim.rphast(obj, c(.001,.002), lower=c(-4,-2), upper=c(2,2)))
system.time(ans[['sym']] <- optim.rphast(obj, tmlol=M[['sym']], c(-1,.4), lower=c(-4,-2), upper=c(2,2)))
objNull <- function(x) obj(c(x, 0))
system.time(ans[['null']] <- optimize(objNull, interval=c(-4,2), maximum=TRUE))

#' ## Simulations of trees conditional on a sampling configuration

getSampleConfig <- function(nm, levs, genTime=as.numeric(ddays(7))){
    splt <- strsplit(nm, split='_')
    tmpf <- function(x) match(x[[1]], levs)
    popIds <- sapply(splt, tmpf)
    tmpf <- function(x) dmy(x[[2]])
    tipDates <- sapply(splt, tmpf)
    sTimes <- max(tipDates) - tipDates
    sGens <- floor(sTimes / genTime)
    tab <- table(sGens, popIds)
    nSamps <- rep(1, length=length(sGens))
    list(ns=nSamps, sg=sGens, sp=popIds)
}

treesim <- function(N=100, nSamples=2, samplingGens=0, samplingPops=1,
                    migProbs=matrix(1, ncol=1,nrow=1), sampleLabels=NULL){
    nPops <- length(N)
    stopifnot(all(samplingPops %in% 1:nPops))
    stopifnot(all(dim(migProbs) == nPops))
    stopifnot(length(samplingGens) == length(samplingPops))
    stopifnot(length(samplingGens) == length(nSamples))
    invMigProbs <- solve(migProbs)
    gen <- min(samplingGens)
    lastSampling <- max(samplingGens)
    totSamples <- sum(nSamples)
    ninternal <- totSamples - 1
    nnode <- totSamples + ninternal
    nedge <- 2*ninternal
    nodeHeights <- numeric(nnode)
    linPops <- integer(nnode)
    nodePops <- numeric(nnode)
    linIds <- 1:nnode
    edge <- matrix(nrow=nedge, ncol=2)
    edge.length <- numeric(nedge)
    nextNode <- nnode
    nextEdge <- 1
    if(!is.null(sampleLabels)){
        stopifnot(length(sampleLabels) == totSamples)
        stopifnot(nSamples == 1)
        tip.label <- character(totSamples)
    }
    isFuture <- rep(FALSE, length.out=nnode)
    nextSample <- 1
    inds <- which(samplingGens == gen)
    for (ind in inds){
        for (i in 1:nSamples[ind]){
            isFuture[nextSample] <- TRUE
            nodeHeights[nextSample] <- gen
            linPops[nextSample] <- samplingPops[ind]
            nodePops[nextSample] <- samplingPops[ind]
            if(!is.null(sampleLabels)){
                tip.label[nextSample] <- sampleLabels[ind]
            }
            nextSample <- nextSample + 1
        }
    }
    isCurrent <- isFuture
    nlin <- sum(isFuture)
    ltt <- list(lineages=list(nlin), jumpTimes=list(gen))
    Nt <- N
    Ntraj <- list(Nt)
    while(nlin > 1 || gen <= lastSampling){
        Ntt <- Nt %*% invMigProbs
        Ntt[Ntt<1] <- 1
        Nt <- Ntt*sum(Nt)/sum(Ntt)
        if (gen %% 10 == 0) {
            Ntraj <- c(Ntraj, list(as.numeric(Nt)))
        }
        gen <- gen + 1
        curPops <- unique(linPops[isCurrent])
        for (pop in curPops) {
            test <- linPops[isCurrent] == pop
            prob <- Nt * migProbs[, pop]
            size <- sum(test)
            if( is.na(size))  browser()
            if( length(prob) != nPops) browser()
            if( any(prob <=0)) browser()
            ancestorPops <- sample(nPops, size=size, replace=TRUE, prob=prob)
            linPops[isCurrent][test] <- ancestorPops
        }
        curPops <- unique(linPops[isCurrent])
        isCurrent <- isFuture
        for (pop in curPops) {
            testPop <- linPops == pop & isCurrent
            nlinPop <- sum(testPop)
            if(nlinPop > 1) {
                ancestors <- sample(ceiling(Nt[pop]), size=nlinPop, replace=TRUE)
                for(i in 1:(nlinPop - 1)){
                    for(j in (i + 1):nlinPop){
                        if (ancestors[i] == ancestors[j]){
                            nodeHeights[nextNode] <- gen
                            lini <- linIds[testPop][i]
                            linj <- linIds[testPop][j]
                            test <- c(lini, linj) %in% edge[, 2]
                            if (test[1] && !test[2]){
                                ind <- which(edge[,2] == lini)
                                parent <- edge[ind, 1]
                                edge[nextEdge, ] <- c(parent, linj)
                                edge.length[nextEdge] <- nodeHeights[parent] - nodeHeights[linj]
                                nextEdge <- nextEdge + 1
                                isFuture[linj] <- FALSE
                            } else if (!test[1] && test[2]){
                                ind <- which(edge[,2] == linj)
                                parent <- edge[ind, 1]
                                edge[nextEdge, ] <- c(parent, lini)
                                edge.length[nextEdge] <- nodeHeights[parent] - nodeHeights[lini]
                                nextEdge <- nextEdge + 1
                                isFuture[lini] <- FALSE
                            } else if (!any(test)){
                                edge[nextEdge, ] <- c(nextNode, lini)
                                edge.length[nextEdge] <- nodeHeights[nextNode] - nodeHeights[lini]
                                nextEdge <- nextEdge + 1
                                edge[nextEdge, ] <- c(nextNode, linj)
                                edge.length[nextEdge] <- nodeHeights[nextNode] - nodeHeights[linj]
                                nextEdge <- nextEdge + 1
                                isFuture[c(lini, linj)] <- FALSE
                                isFuture[nextNode] <- TRUE
                                nodePops[nextNode] <- pop
                                linPops[nextNode] <- pop
                                nextNode <- nextNode - 1
                            }
                        }
                    }
                }
            }
        }
        if (gen %in% samplingGens){
            inds <- which(samplingGens == gen)
            for (ind in inds){
                for (i in 1:nSamples[ind]){
                    isFuture[nextSample] <- TRUE
                    nodeHeights[nextSample] <- gen
                    linPops[nextSample] <- samplingPops[ind]
                    nodePops[nextSample] <- samplingPops[ind]
                    if(!is.null(sampleLabels)) {
                        tip.label[nextSample] <- sampleLabels[ind]
                    }
                    nextSample <- nextSample + 1
                }
            }
        }
        nlinNext <- sum(isFuture)
        if(nlinNext != nlin) {
            ltt$lineages <- c(ltt$lineages, nlinNext)
            ltt$jumpTimes <- c(ltt$jumpTimes, gen)
            nlin <- nlinNext
        }
    }
    if (nextNode != totSamples) {
        diff <- nextNode - totSamples
        ninternal <- ninternal - diff
        nedge <- nedge - diff
        edge <- edge[1:nedge, ]
        edge[edge > totSamples] <- edge[edge > totSamples] - diff
        edge.length <- edge.length[1:nedge]
        start <- totSamples + diff + 1
        end <- start + ninternal - 1
        nodePops <- c(nodePops[1:totSamples], nodePops[start:end])
    }
    if(is.null(sampleLabels)){
        tip.label <- as.character(1:totSamples)
    }
    res <- list(edge=edge,  Nnode=as.integer(ninternal),
                tip.label=tip.label,
                edge.length=edge.length)
    class(res) <- 'phylo'
    ltt <- lapply(ltt, unlist)
    list(phy=res, ltt=ltt, nodePops=nodePops, nodeHeights=nodeHeights, Ntraj=Ntraj)
}

coalStats <- function(ltt){
    jumps <- diff(ltt$lineages)
    holdingTimes <- diff(ltt$jumpTimes)
    k <- ltt$lineages[-length(ltt$lineages)]
    npairs <- k*(k-1)/2
    ci <- list(); cr <- list()
    i <- 1
    while(i < length(jumps)){
        if(jumps[i] < 0) {
            ci <- c(ci, holdingTimes[i])
            cr <- c(cr, npairs[i])
            i <- i + 1
        } else {
            stopifnot(jumps[i] > 0)
            ht <- holdingTimes[c(i, i + 1)]
            np <- npairs[c(i, i + 1)]
            holdingTimes[i + 1] <- sum(ht)
            npairs[i + 1] <- sum(ht*np/sum(ht))
            i <- i + 1
        }
    }
    list(ci=unlist(ci), cr=unlist(cr))
}

get.ltt.phylo <- function(phy){
    nd <- node.depth.edgelength(phy)
    nd <- max(nd) - nd
    nt <- length(phy$tip.label)
    delta <- as.integer(c(rep(1, nt), 1 - table(phy$edge[,1])))
    jumps <- tapply(delta, nd, sum)
    test <- jumps != 0
    lin <- cumsum(jumps[test])
    ltt <- list()
    ltt$lineages <- unname(lin)
    ltt$jumpTimes <- as.numeric(names(lin))
    ltt
}

ran.gen.tree <- function(data=M[['asym']], pars, msal=pedvMSA, levnames=levs,
                         genTime=as.numeric(ddays(7)), branchUnits=as.numeric(ddays(365)),
                         unitsToEvenDist=4, N=70){
    K <- length(levnames)
    gensPerUnit <- branchUnits/genTime
    Z <- data[[1]][[1]]$design.matrix
    Q <- getRateMatrix(Z, pars)
    sampleLabels <- lapply(pedvMSA, names)
    sampCfgs <- lapply(sampleLabels, getSampleConfig, levs=levnames, genTime=genTime)
    migProbs <- expm(Q/gensPerUnit)
    Ni <- N/K
    Ni <- rep(Ni, K) %*% expm(Q*unitsToEvenDist)
    Ni[Ni < 1] <- 1
    Ni <- Ni * N / sum(Ni)
    tmpf <- function(sampCfg, sampleLab){
        treesim(Ni, nSamples=sampCfg$ns, samplingGens=sampCfg$sg,samplingPops=sampCfg$sp, migProbs=migProbs, sampleLabels=sampleLab)
    }
    p <- mapply(tmpf, sampCfgs, sampleLabels, SIMPLIFY=FALSE)
    res <- lapply(data, '[', 1)
    tmpf <- function(x, y) {
        tree <- y$phy
        tree$edge.length <- tree$edge.length / gensPerUnit
        tree <- multi2di(tree)
        x[[1]]$tree <- write.tree(tree)
        x
    }
    res <- mapply(tmpf, x=res, y=p, SIMPLIFY=FALSE)
    res
}

get.param.stat.tree <- function(data, pars) {
    ans <- optim.rphast(obj, params=pars, lower=c(-5,-5), upper=c(2,2), tmlol=data)
    tmpf <- function(x) read.tree(text=x[[1]]$tree)
    trees <- lapply(data, tmpf)
    ltt <- lapply(trees, get.ltt.phylo)
    cs <- lapply(ltt, coalStats)
    tmpf <- function(c){
        y <- c$ci * c$cr
        gini.exp.test(y, simulate=length(y < 30))
    }
    htl <- lapply(cs, tmpf)
    expPvals <- sapply(htl, '[[', 'p.value')
    rej <- expPvals < 0.05
    expStats <- sapply(htl, '[[', 'statistic')
    c(ans$par, rej, expStats)
}

bsParamR <- 1e2
system.time(bsTree <- boot(data=M[['asym']], get.param.stat.tree, R=bsParamR, sim='parametric',
                           ran.gen=ran.gen.tree, mle=ans[['asym']]$par, pars=ans[['asym']]$par,
                           parallel='multicore', ncpus=parallel::detectCores()))

#' ### Bootstrap bias and standard error estimates

bsTree

#' The distributions appears close to normal, without any discontinuities
plot(bsTree, index=1)
plot(bsTree, index=2)
plot(bsTree, index=5)
plot(bsTree, index=6)

#' ## Regularized models and cross-validation

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

dtlmnet <- function(x, y, alpha=1, nlambda=100, lambda.min.ratio=0.01,
                    lambda=NULL, standardize=TRUE, intercept=TRUE, thresh=1e-4,
                    dfmax=nvars + 1, pmax=min(dfmax*2 + 20,nvars), exclude,
                    penalty.factor=rep(1, nvars), lower.limits=-Inf,
                    upper.limits=Inf, maxit=100){
    # written using glmnet as a template
    if (alpha > 1) {
        warning("alpha >1; set to 1")
        alpha <- 1
    }
    if (alpha < 0) {
        warning("alpha<0; set to 0")
        alpha <- 0
    }
    alpha <- as.double(alpha)
    this.call <- match.call()
    nlam <- as.integer(nlambda)
    np <- dim(x)
    if (is.null(np) | (np[2] < 1))
        stop("x should be a matrix with 1 or more columns")
    nrates <- as.integer(np[1])
    nvars <- as.integer(np[2])
    k <- nchar(y$tmlol[[1]][[1]]$alphabet)
    if (k*(k-1) != nrates)
        stop(paste("number of ordered pairs of states for trait (", k*(k-1), ") not consistent with the number of predicted transition rates (", nrates, ")", sep = ""))
    vnames <- colnames(x)
    if (is.null(vnames))
        vnames <- paste("V", seq(nvars), sep = "")
    ne <- as.integer(dfmax)
    nx <- as.integer(pmax)
    if (!missing(exclude)) {
        jd <- match(exclude, seq(nvars), 0)
        if (!all(jd > 0))
            stop("Some excluded variables out of range")
        jd <- as.integer(c(length(jd), jd))
    }
    else jd <- as.integer(0)
    vp <- as.double(penalty.factor)
    if (any(lower.limits > 0)) {
        stop("Lower limits should be non-positive")
    }
    if (any(upper.limits < 0)) {
        stop("Upper limits should be non-negative")
    }
    if (length(lower.limits) < nvars) {
        if (length(lower.limits) == 1)
            lower.limits <- rep(lower.limits, nvars)
        else stop("Require length 1 or nvars lower.limits")
    }
    else lower.limits <- lower.limits[seq(nvars)]
    if (length(upper.limits) < nvars) {
        if (length(upper.limits) == 1)
            upper.limits = rep(upper.limits, nvars)
        else stop("Require length 1 or nvars upper.limits")
    }
    else upper.limits <- upper.limits[seq(nvars)]
    cl = rbind(lower.limits, upper.limits)
    if (any(abs(cl) < thresh)) {
        stop("Cannot enforce limits this close to zero")
    }
    storage.mode(cl) <- "double"
    isd <- as.integer(standardize)
    intr <- as.integer(intercept)
    thresh <- as.double(thresh)
    if (is.null(lambda)) {
        if (lambda.min.ratio >= 1)
            stop("lambda.min.ratio should be less than 1")
        flmin <- as.double(lambda.min.ratio)
        ulam <- double(1)
    }
    else {
        flmin <- as.double(1)
        if (any(lambda < 0))
            stop("lambdas should be non-negative")
        ulam <- as.double(rev(sort(lambda)))
        nlam <- as.integer(length(lambda))
    }
    fit <- dtnet(x, y, alpha, nobs, nvars, jd, vp, cl, ne, nx, nlam, flmin,
                 ulam, thresh, isd, intr, vnames, maxit)
    fit$call <- this.call
    fit$nrates <- nrates
    class(fit) <- c(class(fit), "dtlmnet")
    fit
}

dtnet <- function(x, y, alpha, nobs, nvars, jd, vp, cl, ne, nx, nlam, flmin,
                  ulam, thresh, isd, intr, vnames, maxit, a=0.1, r=0.01,
                  relStart=1, mubar=1, beta=0.9, verbose=TRUE, debug=TRUE,
                  initFactor=10){
    maxit <- as.integer(maxit)
    niter <- 0
    if(!intr) stop("not implemented")
    dim <- nvars + intr
    parInds <- seq(dim)
    I <- diag(nrow=dim)
    nll <- function(w){
        -obj2(w=w, x=x, y=y)
    }
    scaleEst <- getInit(msal=y$msal, tmlol=y$tmlol)
    par <- c(scaleEst, rep(0, nvars))
    gnll <- grad(nll, x=par, method='simple')
    mu <- mubar
    stopifnot(beta>0, beta<1)
    G <- diag(initFactor * abs(gnll), ncol=dim)
    if(flmin<1){
        lstart <- max(abs(gnll))
        loglstart <- log10(lstart) + relStart
        log10LambdaRange <- -log10(flmin)
        loglend <- loglstart - log10LambdaRange
        logstep <- (loglend - loglstart)/nlam
        lambda <- loglstart + 0:(nlam - 1)*logstep
        lambda <- 10^lambda
    } else {
        lambda <- ulam
    }
    res <- list()
    fsg <- function(p, g, l1, l2, h) {
        if(p==0) {
            max(abs(g) - l1, 0)
        } else if (p > 0){
            (-g - l1 - l2*p)/(h + l2)
        } else {
            (-g + l1 - l2*p)/(h + l2)
        }
    }
    for (i in seq_along(lambda)){
        l1penalty <- lambda[i] * c(0, vp) * alpha
        l2penalty <- lambda[i] * c(0, vp) * (1 - alpha)
        penFunc <- function(par) abs(par) %*% l1penalty + (par^2 %*% l2penalty)/2
        k <- 0
        nlp <- nll(par)
        F1 <- nlp + penFunc(par)
        H <- I/(2*mu) + G
        sg <- mapply(fsg, p=par, g=gnll, l1=l1penalty, l2=l2penalty, h=diag(H))
        nsg <- max(abs(sg))
        while (nsg > thresh && k < maxit){
            H <- I/(2*mu) + G
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
                gr <- gnll[j] + Hd[j]
                if (par[j] + d[j] > 0 || (par[j] + d[j] == 0 & -gr > 0)){
                    z <- (-gr - l1penalty[j] - l2penalty[j]*(par[j] + d[j]))/(H[j,j] + l2penalty[j])
                    if (par[j] + d[j] + z < 0){
                        d[j] <- -par[j]
                    } else {
                        d[j] <- d[j] + z
                    }
                } else {
                    z <- (-gr + l1penalty[j] - l2penalty[j]*(par[j] + d[j]))/(H[j,j] + l2penalty[j])
                    if (par[j] + d[j] + z > 0){
                        d[j] <- -par[j]
                    } else {
                        d[j] <- d[j] + z
                    }
                }
                dlist[[nd]] <- d
                fmlist[[nd]] <- d %*% (H/2) %*% d + gnll %*% d + F1 - penFunc(par) + penFunc(par + d)
            }
            if(any(diff(c(F1, unlist(fmlist)))>1e-10)) browser()
            par2 <- par + d
            if (any(par2 > cl[2, ] || any(par2 < cl[1, ]))){
                if(verbose) cat('backtracking: out of bounds', '\n')
                mu <- mu * beta
            } else {
                Fmod <- d %*% (H/2) %*% d + gnll %*% d + F1 - penFunc(par) + penFunc(par2)
                if(Fmod  > F1 && !isTRUE(all.equal(Fmod, F1))) {
                    if (debug) browser()
                    if(verbose) cat('backtracking: poorly solved model', '\n')
                    mu <- mu * beta
                } else {
                    nlp <- nll(par2)
                    F2 <- nlp + penFunc(par2)
                    if (F2 - F1 > r * (Fmod - F1)){
                        if(verbose) cat('backtracking: insufficient decrease', '\n')
                        mu <- mu * beta
                    } else {
                        gnll2 <- grad(nll, x=par2, method='simple')
                        yvec <- gnll2 - gnll
                        s <- d
                        ys <- yvec %*% s
                        if (ys > 0){
                            ## Hessian approximation update via BFGS via 8.19 in Nocedal and Wright
                            sColVec <- matrix(s, ncol=1)
                            yColVec <- matrix(yvec, ncol=1)
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
                        gnll <- gnll2
                        dF <- F2 - F1
                        F1 <- F2
                        H <- I/(2*mu) + G
                        sg <- mapply(fsg, p=par, g=gnll, l1=l1penalty, l2=l2penalty, h=diag(H))
                        nsg <- max(abs(sg))
                        if (debug && mu < 1e-8 && nsg > thresh) browser()
                        if (verbose) {
                            cat('lambda: ', lambda[i], '\n')
                            cat('k: ', k, '\n')
                            cat('F: ', as.numeric(F1), '\n')
                            cat('par: ', signif(par, 3), '\n')
                            cat('grad: ', signif(gnll, 3), '\n')
                            cat('nsg: ', signif(nsg, 3), '\n')
                            cat('\n')
                        }
                    }
                }
            }
        }
        convergence <- ifelse(k == maxit, 'no', 'yes')
        res[[i]] <- list(par=par, nll=nlp, k=k, lambda=lambda[i], convergence=convergence, mu=mu, nsg=nsg, sg=sg)
        mu <- mubar
    }
    path <- sapply(res, '[[', 'par')
    ret <- list(a0=path[1,])
    beta <- t(path[-1, ])
    colnames(beta) <- paste("s", seq(ncol(beta)) - 1, sep = "")
    rownames(beta) <- vnames
    ret$beta <- beta
    ret$lambda <- sapply(res, '[[', 'lambda')
    ret$nll <- sapply(res, '[[', 'nll')
    ret$df <- colSums(beta > 0)
    ret$dim <- dim(x)
    ret$niterations <- sum(sapply(res, '[[', 'k'))
    ret$jerr <- paste('convergence: ', sapply(res, '[[', 'convergence'))
    ret
}

plotCoef <- function (beta, norm, lambda, df, dev, label = FALSE, xvar = c("norm",
    "lambda", "dev"), xlab = iname, ylab = "Coefficients", ...) {
    which <- which(rowSums(abs(beta)) > 0)
    nwhich <- length(which)
    switch(nwhich + 1, `0` = {
        warning("No plot produced since all coefficients zero")
        return()
    }, `1` = warning("1 or less nonzero coefficients; glmnet plot is not meaningful"))
    beta = as.matrix(beta[which, , drop = FALSE])
    xvar = match.arg(xvar)
    switch(xvar, norm = {
        index = if (missing(norm)) apply(abs(beta), 2, sum) else norm
        iname = "L1 Norm"
        approx.f = 1
    }, lambda = {
        index = log(lambda)
        iname = "Log Lambda"
        approx.f = 0
    }, dev = {
        index = dev
        iname = "Fraction Deviance Explained"
        approx.f = 1
    })
    dotlist = list(...)
    type = dotlist$type
    if (is.null(type))
        matplot(index, t(beta), lty = 1, xlab = xlab, ylab = ylab,
            type = "l", ...)
    else matplot(index, t(beta), lty = 1, xlab = xlab, ylab = ylab,
        ...)
    atdf = pretty(index)
    prettydf = approx(x = index, y = df, xout = atdf, rule = 2,
        method = "constant", f = approx.f)$y
    axis(3, at = atdf, labels = prettydf, tcl = NA)
    if (label) {
        nnz = length(which)
        xpos = max(index)
        pos = 4
        if (xvar == "lambda") {
            xpos = min(index)
            pos = 2
        }
        xpos = rep(xpos, nnz)
        ypos = beta[, ncol(beta)]
        text(xpos, ypos, paste(which), cex = 0.5, pos = pos)
    }
}

plot.dtlmnet <- function (x, xvar = c("norm", "lambda", "dev"),
                          label = FALSE, ...){
    xvar = match.arg(xvar)
    plotCoef(x$beta, lambda = x$lambda, df = x$df, dev = x$dev.ratio,
             label = label, xvar = xvar, ...)
}


x <- M[[1]][[1]][[1]]$design.matrix
y <- list(tmlol=M[['big']], msal=pedvMSA)
dfit <- dtlmnet(x=x, y=y)
plot(dfit, xvar='l')



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

save.image('mle.RData')
