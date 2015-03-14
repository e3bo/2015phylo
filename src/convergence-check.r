library(coda)

burnin <- 0.1
runstems <- c('pedv', 'pedv-run2')

files <- paste(runstems, '.log', sep='')
runs <- lapply(files, read.table, header=TRUE, sep="\t")

tmpf <- function(x) {
    diff(x$state[1:2])
}
thin <- sapply(runs, tmpf)
stopifnot(all.equal(thin, rep(thin[1], length.out=length(thin))))
thin <- thin[1]

tmpf <- function(x) {
    x[, 'state'] <- NULL
    x
}
runs <- lapply(runs, tmpf)

mobs <- lapply(runs, mcmc, thin=thin)
ml <- mcmc.list(mobs)
wstart <- burnin*end(ml)

wml <- window(ml, start=wstart)
ess <- lapply(wml, effectiveSize)

gd <- gelman.diag(wml, multivariate=FALSE)
