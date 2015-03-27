#!/usr/bin/Rscript

library(ape)
library(coda)
library(ggplot2)
library(lubridate)
library(reshape2)

burnin <- 0.5
runstems <- paste0('beast/run', 1:4, '/pedv')

files <- paste(runstems, '.log', sep='')
runs <- lapply(files, read.table, header=TRUE, sep="\t")
prior <- read.table('beast/priors/pedv-priors.log', header=TRUE, sep='\t')
runsPrior <- c(runs, list(prior))

tmpf <- function(x) {
    diff(x$state[1:2])
}
thin <- sapply(runsPrior, tmpf)
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
runsp <- lapply(runsPrior, tmpf)

mobs <- lapply(runsp, mcmc, thin=thin)
ml <- mcmc.list(mobs[seq_along(runs)])
wstart <- burnin*end(ml)
wml <- window(ml, start=wstart)

i <- length(runs) + 1
pr <- mobs[[i]]
wstart <- burnin*end(pr)
wpr <- window(pr, start=wstart)

plotVars <- c("nonsIndel.aligned.fasta.gb.treeModel.rootHeight",
              "sIndel.aligned.fasta.gb.treeModel.rootHeight")
mml <- as.matrix(wml[, plotVars])
dpst <- data.frame(mml, Function='posterior')

mpr <- as.matrix(wpr[, plotVars])
dpr <- data.frame(mpr, Function='prior')

D <- rbind(dpst, dpr)
colnames(D) <- gsub('.aligned.fasta.gb.treeModel.rootHeight', '', colnames(D))
colnames(D) <- gsub('sIndel', 'S INDEL', colnames(D))
colnames(D) <- gsub('nonS', 'non-S', colnames(D))
M <- melt(D, id.vars='Function')

theme_set(theme_classic())

files <- c('nonsIndel-aligned.fasta-gb', 'sIndel-aligned.fasta-gb')
tmpf <- function(x){
    x <- read.dna(x, format='fasta')
    rownames(x)
}
snames <- unname(unlist(sapply(files, tmpf)))
snames <- strsplit(snames, '_')
sdates <- sapply(snames, '[[', 2)
sdates <- dmy(sdates)
tzero <- max(sdates)
xlabel <- paste("\nYears before", format(max(sdates), "%B %Y"))

firstOutbreak <- dmy("15-Apr-2013")
firstOutbreakYearsBeforeTzero <- as.duration(tzero - firstOutbreak)/dyears(1)

colScale <- c("prior"="#4477AA", "posterior"="#CC6677")

g <- ggplot(data=M, aes(x=value, color=Function))
g <- g + theme(legend.position='top')
g <- g + geom_vline(xintercept=firstOutbreakYearsBeforeTzero, col='grey')
g <- g + geom_density()
g <- g + scale_colour_manual(values=colScale)
g <- g + xlim(c(0,5))
g <- g + xlab(xlabel) + ylab('Density\n')
g <- g + facet_grid(variable~.)
ggsave('posteriors-and-priors.pdf', width=3.5)
