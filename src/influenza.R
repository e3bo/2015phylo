
library(ape)
raxmlbin <- "/usr/bin/raxmlHPC"
alignbin <- "../data/swine-influenza-alignments-dnabin.rds"

stopifnot(file.exists(raxmlbin))
stopifnot(file.exists(alignbin))


align <- readRDS(alignbin)

tr <- ips::raxml(align$H1N1, m = "GTRGAMMA", p = 12345, N = 3, f = "a", exec = raxmlbin)
rate_ests <- get_raxml_ests(tr = tr)
bt <- tr$bestTree

