library(ape)
library(coda)

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
(gd <- gelman.diag(wml, multivariate=FALSE))
gelman.diag(wml, transform=TRUE, multivariate=FALSE)

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
my.gelman.plot(wml, ylim=c(1, 1.1))

save.image('convergence-check.RData')


tfiles <- paste(runstems, '.trees', sep='')
truns <- lapply(tfiles, read.nexus)

ntrees <- sapply(truns, length)
stopifnot(isTRUE(all.elements.equal(ntrees)))
ntrees <- ntrees[1]

brn <- seq_len(ceiling(ntrees*burnin))
tmpf <- function(x) {
    x[-brn]
}
wtrns <- lapply(truns, tmpf)

prt <- lapply(wtrns, prop.part)
