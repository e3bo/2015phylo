#!/usr/bin/Rscript

load("mers-c1.RData")

x2 <- diag(4)[, -1]
tree_timel <- lapply(tree_info, "[[", "tree_time")
KSAtips <- grep("Riyadh|Taif|Jeddah|KSA", tree_timel[[1]]$tip.label)
KSAtips <- tree_timel[[1]]$tip.label[KSAtips]

pm2 <- function (x, w){
  psampled <- c(.01, .1)
  n <- 2
  ret <- list(frequency = list(c(1, 0)))
  ret$m <- rep(exp(w[1]), n) * (1 - psampled)
  ret$psi <- rep(exp(w[1]), n) * psampled
  scale <- exp(w[2])
  effects <- w[-c(1,2)]
  stopifnot(nrow(x) == n^2)
  stopifnot(ncol(x) == length(effects))
  eta <- exp(x %*% effects)
  eta <- eta/mean(eta) * scale
  rate_matrix <- matrix(eta, nrow = n, ncol = n)
  ret$l <- rate_matrix
  ret$survival <- FALSE
  ret
}                   
                                
                                
init1 <- c(2, 1, rep(0,3))
pf1 <- c(0, 0, rep(1, 3))

pars <- pm2(x = x2, w = init1)

tree_timel <- lapply(tree_info, "[[", "tree_time")

tree_timelf <- penaltree:::filter_y(tree_timel, keepers = KSAtips)

ncpu <- min(parallel::detectCores() - 1, 30)
sp <- penaltree::stabpath_gpnet(x = x2, y = tree_timelf[1],
               calc_convex_nll = penaltree::calc_bd_lm_nll,
               param_map = pm2, nlambda = 10, lambda.min.ratio = 0.1,
               make_log = TRUE, penalty.factor = pf1,
               thresh = 1e-3, winit = init1, alpha = 1,
               steps=30, mc.cores = ncpu)
save.image("mers-c2.RData")
