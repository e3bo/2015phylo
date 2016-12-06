#!/usr/bin/Rscript

library(c060)
load("influenza-c2.RData")

spstats <- plot(sp)
xstable <- x2[, 4L - 2L, drop=FALSE]
stopifnot("V4" %in% names(spstats$stable))

pf2 <- c(0, 0, 100, 100, 1)
init2 <- init1[c(seq(1, 4), 4L + 4L)]
selected_fit <- penaltree::get_gpnet(x = xstable, y = tree_timel[1],
               calc_convex_nll = penaltree::calc_bd_lm_nll,
               param_map = pm1, nlambda = 10, lambda.min.ratio = 0.01,
               make_log = TRUE, penalty.factor = pf2,
               thresh = 1e-3, winit = init2, alpha = 1)
save.image("influenza-c3.RData")
