library(ape)
library(coda)

my.gelman.preplot <- function (x, bin.width = bin.width, max.bins = max.bins, confidence = confidence, 
    transform = transform, autoburnin = autoburnin) 
{
    ## based on code in coda:::gelman.preplot
    ## this version avoids calculating the multivariate PSRF, which lead to errors for some bins
    x <- as.mcmc.list(x)
    if (niter(x) <= 50) 
        stop("Less than 50 iterations in chain")
    nbin <- min(floor((niter(x) - 50)*thin(x)/bin.width), max.bins)
    binw <- floor((niter(x) - 50)/nbin)
    last.iter <- c(seq(from = start(x) + 50 * thin(x), by = binw * 
        thin(x), length = nbin), end(x))
    shrink <- array(dim = c(nbin + 1, nvar(x), 2))
    dimnames(shrink) <- list(last.iter, varnames(x), c("median", 
        paste(50 * (confidence + 1), "%", sep = "")))
    for (i in 1:(nbin + 1)) {
        shrink[i, , ] <- gelman.diag(window(x, end = last.iter[i]), 
            confidence = confidence, transform = transform, autoburnin = autoburnin, multivariate=F)$psrf
    }
    all.na <- apply(is.na(shrink[, , 1, drop = FALSE]), 2, all)
    if (any(all.na)) {
        cat("\n******* Error: *******\n")
        cat("Cannot compute Gelman & Rubin's diagnostic for any chain \n")
        cat("segments for variables", varnames(x)[all.na], "\n")
        cat("This indicates convergence failure\n")
    }
    return(list(shrink = shrink, last.iter = last.iter))
}

my.gelman.plot <- function (x, bin.width = 10, max.bins = 50, confidence = 0.95, 
    transform = FALSE, autoburnin = TRUE, auto.layout = TRUE, 
    ask, col = 1:2, lty = 1:2, xlab = "last iteration in chain", 
    ylab = "shrink factor", type = "l", ...) 
{
    ## same as coda::gelman.plot but uses my.gelman.preplot
    if (missing(ask)) {
        ask <- if (is.R()) {
            dev.interactive()
        }
        else {
            interactive()
        }
    }
    x <- as.mcmc.list(x)
    oldpar <- NULL
    on.exit(par(oldpar))
    if (auto.layout) 
        oldpar <- par(mfrow = coda:::set.mfrow(Nchains = nchain(x), 
            Nparms = nvar(x)))
    y <- my.gelman.preplot(x, bin.width = bin.width, max.bins = max.bins, 
        confidence = confidence, transform = transform, autoburnin = autoburnin)
    all.na <- apply(is.na(y$shrink[, , 1, drop = FALSE]), 2, 
        all)
    if (!any(all.na)) 
        for (j in 1:nvar(x)) {
            matplot(y$last.iter, y$shrink[, j, ], col = col, 
                lty = lty, xlab = xlab, ylab = ylab, type = type, 
                ...)
            abline(h = 1)
            ymax <- max(c(1, y$shrink[, j, ]), na.rm = TRUE)
            leg <- dimnames(y$shrink)[[3]]
            xmax <- max(y$last.iter)
            legend(xmax, ymax, legend = leg, lty = lty, bty = "n", 
                col = col, xjust = 1, yjust = 1)
            title(main = varnames(x)[j])
            if (j == 1) 
                oldpar <- c(oldpar, par(ask = ask))
        }
    return(invisible(y))
}

burnin <- 0.5
runstems <- paste0('simplest/run', 1:4, '/simplest')

files <- paste(runstems, '.log', sep='')
runs <- lapply(files, read.table, header=TRUE, sep="\t")

tmpf <- function(x) {
    diff(x$state[1:2])
}
thin <- sapply(runs, tmpf)
all.elements.equal <- function(x) all.equal(x, rep(x[1], length.out=length(x)))
stopifnot(isTRUE(all.elements.equal(thin)))
thin <- thin[1]

tmpf <- function(x) {
    x[, 'state'] <- NULL
#    x[, 'clock.rate'] <- NULL
#    x[, 'posterior'] <- NULL
#    x[, 'coalescent'] <- NULL
#    x[, 'treeLikelihood'] <- NULL
#    x[, 'prior'] <- NULL
    x
}
runsp <- lapply(runs, tmpf)

mobs <- lapply(runsp, mcmc, thin=thin)
ml <- mcmc.list(mobs)
wstart <- burnin*end(ml)
wml <- window(ml, start=wstart)

plot(wml)

(essl <- lapply(wml, effectiveSize))
(ess <- effectiveSize(wml))

heidel.diag(wml)

(gd <- gelman.diag(wml, multivariate=FALSE))
gelman.diag(wml, transform=TRUE, multivariate=FALSE)

#' ## Gelman plots
my.gelman.plot(wml, ylim=c(1, 1.1))

#' ## Geweke plots
geweke.plot(wml)

#' ##Now for the trees

phymcmc <- function(x, start=1, end=numeric(0), thin=1){
    ## extension of mcmc to multiPhylo objects, based on code in coda 
    if (!inherits(x, 'multiPhylo')) {
        stop("Not an object of class multiPhylo")
    } else {
        niter <- length(x)
        nms <- names(x)[1:2]
        s12 <- as.integer(sapply(strsplit(nms, split='_'), '[[', 2))
        if(missing(thin))
            thin <- diff(s12)
        if(missing(start))
            start=s12[1]
        if(missing(end))
            end <- start + thin * (niter - 1)
        nobs <- floor((end - start)/thin + 1)
        if (niter != nobs)
            stop("Start, end and thin incompatible with data")
        attr(x, 'mcpar') <- c(start, end, thin)
        attr(x, 'class') <- c('phymcmc', 'multiPhylo')
        x
    }
}

start.phymcmc <- function (x, ...)
    {
        mcpar(x)[1]
    }

end.phymcmc <- function (x, ...)
    {
        mcpar(x)[2]
    }

frequency.phymcmc <- function (x, ...)
    {
        1/thin(x)
    }

thin.phymcmc <- function (x, ...)
    {
        mcpar(x)[3]
    }

time.phymcmc <- function (x, ...)
    {
        ts(seq(from = start(x), to = end(x), by = thin(x)), start = start(x),
           end = end(x), deltat = thin(x))
    }

window.phymcmc <- function (x, start, end, thin, ...) 
{
    ts.eps <- getOption("ts.eps")
    xmcpar <- mcpar(x)
    xstart <- xmcpar[1]
    xend <- xmcpar[2]
    xthin <- xmcpar[3]
    if (missing(thin)) 
        thin <- xthin
    else if (thin%%xthin != 0) {
        thin <- xthin
        warning("Thin value not changed")
    }
    xtime <- as.vector(time(x))
    if (missing(start)) 
        start <- xstart
    else if (length(start) != 1) 
        stop("bad value for start")
    else if (start < xstart) {
        start <- xstart
        warning("start value not changed")
    }
    if (missing(end)) 
        end <- xend
    else if (length(end) != 1) 
        stop("bad value for end")
    else if (end > xend) {
        end <- xend
        warning("end value not changed")
    }
    if (start > end) 
        stop("start cannot be after end")
    if (all(abs(xtime - start) > abs(start) * ts.eps)) {
        start <- xtime[(xtime > start) & ((start + xthin) > xtime)]
    }
    if (all(abs(end - xtime) > abs(end) * ts.eps)) {
        end <- xtime[(xtime < end) & ((end - xthin) < xtime)]
    }
    use <- 1:length(x)
    use <- use[use >= trunc((start - xstart)/xthin + 1.5) & use <= 
        trunc((end - xstart)/xthin + 1.5) & (use - trunc((start - 
        xstart)/xthin + 1.5))%%(thin%/%xthin) == 0]
    y <- x[use]
    return(phymcmc(y, start = start, end = end, thin = thin))
}

print.phymcmc <- function (x, ...)
    {
        x.orig <- x
        cat("Markov Chain Monte Carlo (MCMC) output:\nStart =", start(x),
            "\nEnd =", end(x), "\nThinning interval =", thin(x), "\n")
        attr(x, "mcpar") <- NULL
        attr(x, "class") <- 'multiPhylo'
        NextMethod("print", ...)
        invisible(x.orig)
    }

is.phymcmc <- function(x) inherits(x, 'phymcmc')
    
phymcmc.list <- function (...)
    {
        x <- list(...)
        if (length(x) == 1 && is.list(x[[1]]))
            x <- x[[1]]
        if (!all(unlist(lapply(x, is.phymcmc))))
            stop("Arguments must be phymcmc objects")
        nargs <- length(x)
        if (nargs >= 2) {
            xmcpar <- lapply(x, mcpar)
            if (!all(unlist(lapply(xmcpar, "==", xmcpar[[1]]))))
                stop("Different start, end or thin values in each chain")
            xtlabels <- lapply(x, attr, 'TipLabel')
            if (!all(unlist(lapply(xtlabels, "==", xtlabels[[1]]))))
                stop("Different taxa among chains")            
        }
        class(x) <- "phymcmc.list"
        return(x)
    }

start.phymcmc.list <- function (x, ...)
    {
        start(x[[1]])
    }

end.phymcmc.list <- function (x, ...)
    {
        end(x[[1]])
    }

thin.phymcmc.list <- function (x, ...)
    {
        thin(x[[1]])
    }

window.phymcmc.list <- function (x, ...)
    {
        structure(lapply(x, window.phymcmc, ...), class = "phymcmc.list")
    }

getSDSF <- function(x, minCladeFz=0.1){
    niter <- length(x[[1]])
    prt <- lapply(x, prop.part)
    labs <- lapply(prt, attr, 'labels')
    allSame <- function(v) all(sapply( v[-1], FUN=function(z) {identical(z, v[[1]])}))
    stopifnot(allSame(labs))
    tmpf <- function(x) {
        count <- attr(x, 'number')
        tmpf <- function(x){
            stopifnot(!is.unsorted(x, strictly=TRUE))
            paste(x, collapse=',')
        }
        clade <- sapply(x, tmpf)
        data.frame(clade=clade, count=count, stringsAsFactors=FALSE)
    }
    foo <- lapply(prt, tmpf)
    tmpf <- function(x, y) {
        n <- ncol(x)
        suff <- as.character(c(n-1, n))
        merge(x, y, by='clade', all=TRUE, suffixes=suff)
    }
    mg <- Reduce(f=tmpf, x=foo)
    cols <- grep('count', colnames(mg))
    counts <- mg[, cols]
    rownames(counts) <- mg$clade
    counts[is.na(counts)] <- 0
    fz <- counts/niter
    test <- apply(fz, 1, function(x) max(x) > minCladeFz)
    fz <- fz[test,]
    sdcf <- apply(fz, 1, sd)
    sdcf
}

asdsf.preplot <- function(x, nbin){
    niter <- length(x[[1]])
    if(niter<=50){
        stop("Less than 50 iterations in chain")
    }
    binw <- floor((niter - 50)/nbin)
    last.iter <- c(seq(from=start(x) + 50*thin(x), by=binw * thin(x) , length=nbin), end(x))
    tmpf <- function(y) {
        z <- window(x, end=y)
        getSDSF(z)
    }
    sdsf <- lapply(last.iter, tmpf)
    return(list(sdsfl=sdsf, last.iter=last.iter))
}

asdsf.plot <- function(...){
    ans <- asdsf.preplot(...)
    means <- sapply(ans$sdsfl, mean)
    plot(c(ans$last.iter[1], ans$last.iter), c(0, means), type='n',
         xlab='Last iteration', ylab='ASDCF')
    lines(ans$last.iter, means)
    abline(h=0, col='grey')
    abline(h=0.01, col='blue', lty=2)
    return(invisible(ans))
}

## Now to actually make use of all those functions

tfiles <- paste(runstems, '.trees', sep='')
truns <- lapply(tfiles, read.nexus)
tmobs <- lapply(truns, phymcmc)
tml <- phymcmc.list(tmobs)

wstart <- burnin*end(tml)
wtml <- window(tml, start=wstart)

nbin <- 8
system.time(ans <- asdsf.plot(wtml, nbin=nbin))
hist(log10(ans$sdsfl[[nbin]]))

save.image('convergence-check.RData')
