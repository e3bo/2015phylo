set_branchlengths <- function(tree, nodeheights, tipheights){
    tree <- reorder(tree, order="postorder")
    edge <- tree$edge
    ntips <- length(tipheights)
    nedge <- nrow(edge)
    stopifnot(nedge + 1 == length(nodeheights) + ntips)
    stopifnot(ntips == length(tree$tip.label))
    nh <- numeric(nedge + 1)
    idx_tip <- which(edge[, 2] <= ntips)
    idx_internal <- c(which(edge[, 2] > ntips), nedge + 1)
    stopifnot(length(setdiff(names(tipheights), tree$tip.label)) == 0)
    tipnames <- tree$tip.label[edge[idx_tip, 2]]
    nh[idx_tip] <- tipheights[tipnames]
    nh[idx_internal] <- nodeheights
    ipheight <- match(edge[, 1], edge[,2])
    ipheight[is.na(ipheight)] <- nrow(edge) + 1
    ## Ensure ordering of nodeheights is valid
    penalty <-0
    for (i in seq(1, nedge)){
        blen <- nh[ipheight[i]] - nh[i]
        if (blen <= 0){
            penalty <- penalty + blen^2
            nh[ipheight[i]] <- nh[i] + .Machine$double.eps
        }
    }
    ## Set branch lengths to differences in heights
    for (i in seq(1, nedge)){
        tree$edge.length[i] <- nh[ipheight[i]] - nh[i]
    }
    list(tree=tree, penalty=penalty)
}

get_time_tree_internal_nodeheights <- function(subs_tree, subs_per_time, tipheights) {
    time_tree <- subs_tree
    time_tree$edge.length <- subs_tree$edge.length / subs_per_time
    est_nh <- get_nodeheights(time_tree)
    delta <- mean(tipheights) - mean(est_nh$tipheights)
    nodeheights <- est_nh$nodeheights + delta
    ifelse(nodeheights < 0, .Machine$double.eps, nodeheights)
}

pml_wrapper <- function(tree, nodeheights, tipheights, phydata, ...){
    foo <- set_branchlengths(tree, nodeheights, tipheights)
    print(nodeheights)
    nll <- -phangorn::pml(foo$tree, data=phydata, ...)$logLik
    print(paste("nll = ", nll))
    print(paste("penalty =", foo$penalty))
    nll + foo$penalty
}

get_hky_Q <- function(kappa=4, pi=c(A=.2, C=.25, G=.3, T=.25)){
    stopifnot(isTRUE(all.equal(sum(pi), 1)))
    pi <- pi[c("A", "C", "G", "T")]
    Q <- matrix(rep(pi, each=4), ncol=4)
    Q[1,3] <- Q[1,3] * kappa
    Q[3,1] <- Q[3,1] * kappa
    Q[2,4] <- Q[2,4] * kappa
    Q[4,2] <- Q[4,2] * kappa
    dimnames(Q) <- list(to=names(pi), from=names(pi))
    diag(Q) <- 0
    diag(Q) <- -rowSums(Q)
    scale <- -sum(diag(Q) * pi)
    Q / scale
}


get_tree_subs <- function(tree_time, node_times, tip_times){
    tree_subs <- set_branchlengths(tree_time, node_times, tip_times)
    tree_subs$tree$edge.length <- tree_subs$tree$edge.length * subs_per_time
    tree_subs
}

lmsa_wrapper <- function(tree_time, node_times, tip_times, msa, subs_per_time,
                         subs_model, alpha, nrates, subs_pars, pi){
    tree_subs <- get_tree_subs(tree_time, node_times, tip_times)
    rate.consts <- phangorn::discrete.gamma(alpha, nrates)
    rate.weights <- rep(1 / nrates, nrates)
    tree_char <- ape::write.tree(tree_subs$tree)
    tmod <- rphast::tm(tree_char, subs_model, backgd = pi, nratecats = nrates,
                       rate.consts = rate.consts, rate.weights = rate.weights)
    tmod <- rphast::set.rate.matrix.tm(tmod, params=subs_pars, scale=TRUE)
    rphast::likelihood.msa(x=msa, tm=tmod) - tree_subs$penalty
}

get_nodeheights <- function(tree, dist_from_root = FALSE){
    tpo <- reorder(tree, "postorder")
    edge <- tpo$edge
    ntips <- length(tpo$tip.label)
    nedge <- nrow(edge)
    nh <- numeric(nedge + 1)
    idx_root <- nedge + 1
    idx_tip <- which(edge[, 2] <= ntips)
    idx_internal <- c(which(edge[, 2] > ntips), idx_root)
    nh[idx_root] <- 0
    ipheight <- match(edge[, 1], edge[,2])
    ipheight[is.na(ipheight)] <- idx_root
    for (i in seq(nedge, 1)){
        nh[i] <- nh[ipheight[i]] - tpo$edge.length[i]
    }
    if (!dist_from_root) {
        nh <- nh - min(nh)
    } else {
        nh <- -nh
    }
    tiph <- nh[idx_tip]
    ind <- edge[idx_tip, 2]
    names(tiph) <- tpo$tip.label[ind]
    list(nodeheights=nh[idx_internal], tipheights=tiph, nh=nh)
}

eval_temporal_signal <- function(phy, tipheights, show_plots=FALSE){
    root_tip_dists <- get_nodeheights(phy, dist_from_root=TRUE)$tipheights
    nms <- names(root_tip_dists)
    tip_times <- -tipheights[nms]
    m <- lm(root_tip_dists ~ tip_times)
    if(show_plots){
        plot(root_tip_dists ~ tip_times)
        abline(m)
        plot(m)
    }
    if (coef(m)[1] < 0){
        print("Slope is not positive, no strong signal")
    }
    ret <- list()
    ret$subs_per_time <- coef(m)[2]
    ret$tmrca <- -coef(m)[1] / coef(m)[2]
    ret$date_range <- range(tipheights)
    sstot <- sum((root_tip_dists - mean(root_tip_dists))^2)
    ssres <- sum(residuals(m)^2)
    ret$coef_det <- 1 - ssres / sstot
    ret$rmes <- ssres / (length(root_tip_dists) - 2)
    ret$mod <- m
    ret
}

get_raxml_ests <- function(tr){
    ret <- list()
    get_freq <- function(char){
        pattern <- paste0("freq pi\\(", char, "\\): [0-9.]*$")
        iline <- grep(pattern, tr$info)
        line <- tr$info[iline]
        as.numeric(strsplit(line, ": ")[[1]][2])
    }
    alphabet <- c("A", "C", "G", "T")
    ret$bf <- sapply(alphabet, get_freq)
    get_alpha <- function(){
        pattern <- paste0("^alpha: [0-9.]*$")
        iline <- grep(pattern, tr$info)
        line <- tr$info[iline]
        as.numeric(strsplit(line, ": ")[[1]][2])
    }
    ret$alpha <- get_alpha()
    get_gtr_pars <- function(alphabet){
        ret <- list()
        ind <- 1
        for (i in seq_along(alphabet)){
            for (s2 in alphabet[-seq(1, i)]) {
                s1 <- alphabet[i]
                name <- paste0(s1, "_", s2)
                pattern <- paste0("^rate ", s1, " <-> ", s2, ": [0-9.]*$")
                iline <- grep(pattern, tr$info)
                if (length(iline) == 0){
                    pattern <- paste0("^rate ", s2, " <-> ", s1, ": [0-9.]*$")
                    iline <- grep(pattern, tr$info)
                }
                line <- tr$info[iline]        
                ret[[name]] <- as.numeric(strsplit(line, ": ")[[1]][2])
            }
        }
        unlist(ret)
    }
    ret$gtr_pars <- get_gtr_pars(alphabet)
    ret
}

#' Generate parameter map for phylogeny nll
#'
#' @export
gen_param_map_phylo <- function(tree, tip_times){
    function(w){
        ret <- list(tree=tree, tip_times=tip_times)
        ret$subs_per_time <- w[1]
        ret$pi <- w[seq(2, 4)]
        ret$pi <- c(ret$pi, 1) / (1 + sum(ret$pi))
        names(ret$pi) <- c("A", "C", "G", "T")
        ret$subs_pars <- c(w[seq(5, 9)], 1)
        ret$alpha <- w[10]
        ret$node_times <- w[seq(11, length(w))]
        ret
    }
}

#' Calculate the negative log likelihood of a time scaled phylogeny
#'
#' @export
calc_phylo_nll <- function(w, x, y, param_map){
    pars <- param_map(w)
    try(nll <- -lmsa_wrapper(pars$tree, node_times = pars$node_times, tip_times = pars$tip_times,
                 msa = y, subs_per_time = pars$subs_per_time, alpha = pars$alpha,
                 subs_model = "REV", nrates = 4L, subs_pars = pars$subs_pars,
                             pi = pars$pi))
    if (inherits(nll, "try-error")){
        browser()
    } else {
        nll
    }
}

#' Generate parameter map for phylogeny nll
#'
#' @export
gen_param_map_phylo_bd <- function(tree, tip_times){
    function(x, w){
        ncoefs <- ncol(x) + 1
        wph <- w[seq(1, length(w) - ncoefs)]
        wbd <- w[seq(length(w) - ncoefs + 1, length(w))]
        ret <- list(tree=tree, tip_times=tip_times)
        ret$subs_per_time <- wph[1]
        ret$pi <- wph[seq(2, 4)]
        ret$pi <- c(ret$pi, 1) / (1 + sum(ret$pi))
        names(ret$pi) <- c("A", "C", "G", "T")
        ret$subs_pars <- c(wph[seq(5, 9)], 1)
        ret$alpha <- wph[10]
        ret$node_times <- wph[seq(11, length(wph))]

        scale <- exp(wbd[1])
        effects <- wbd[-1]
        eta <- exp(x %*% effects)
        eta <- eta / mean(eta) * scale
        n <- sqrt(nrow(x))
        rate_matrix <- matrix(eta, nrow=n, ncol=n)
        ret$l <- rate_matrix
        ret$m <- rep(1, n) / 2
        ret$psi <- rep(1, n) / 2
        ret$survival <- FALSE
        ret$frequency <- c(1, rep(0, n - 2))
        ret
    }
}

#' Calculate the negative log likelihood of a time scaled phylogeny
#' assuming it was created by a birth-death process
#'
#' @export
calc_phylo_nll_bd <- function(w, x, y, param_map){
    pars <- param_map(x, w)
    nllphy <- -lmsa_wrapper(pars$tree, node_times = pars$node_times, tip_times = pars$tip_times,
                 msa = y, subs_per_time = pars$subs_per_time, alpha = pars$alpha,
                 subs_model = "REV", nrates = 4L, subs_pars = pars$subs_pars,
                         pi = pars$pi)
    phylo <- get_tree_subs(pars$tree, pars$node_times, pars$tip_times)$tree
    phylo <- TreePar::addroot(phylo, phylo$root.edge)
    nllbd <- calc_bd_nll(l=pars$l, m=pars$m, psi=pars$psi, freq=pars$freq, phylo=phylo,
                         survival=pars$survival)
    nllphy + nllbd
}
