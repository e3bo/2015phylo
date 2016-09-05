
library(ape)
trees <- read.nexus("pedv.nonsIndel-aligned.fasta-gb.trees")
tree <- trees[[10001]]
checkValidPhylo(tree)

datepos <- regexpr("([0-9]{2}\\-[A-Z][a-z]{2}\\-[0-9]{4})", tree$tip.label)
tmpf <- function(x, y) substring(x, y, y+10)
dates <- mapply(tmpf, tree$tip.label, datepos)
dates <- as.Date(dates, "%d-%b-%Y")

library(lubridate)
decdates <- year(dates) + (month(dates) - 1) / 12 + day(dates) / 365
tree$tip.label <- paste(tree$tip.label, format(decdates, digits=6, nsmall=2), sep="_")

fn <- "tipdate.in"
cat("1\n", file=fn, append=TRUE)
write.tree(tree, file=fn, append=TRUE)
