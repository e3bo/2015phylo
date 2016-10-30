
set_branchlengths <- function(tree, nodeheights, tipheights){
    tree <- reorder(tree, order="postorder")
    edge <- tree$edge
    ntips <- length(tipheights)
    nedge <- nrow(edge)
    stopifnot(nedge + 1== length(nodeheights) + ntips)
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

lmsa_wrapper <- function(tree_time, node_times, tip_times, msa, subs_per_time,
                         subs_model, subs_pars){
    tree_subs <- set_branchlengths(tree_time, node_times, tip_times)
    tree_subs$tree$edge.length <- tree_subs$tree$edge.length * subs_per_time
    tree_char <- ape::write.tree(tree_subs$tree)
    if (subs_model == "HKY85") {
        Q <- get_hky_Q(subs_pars$kappa, subs_pars$pi)
        tmod <- rphast::tm(tree_char, subs_model, rate.matrix=Q, backgd = subs_pars$pi)
    } else {
        stop("subs_model not implemented")
    }
    rphast::likelihood.msa(x=msa, tm=tmod) - tree_subs$penalty
}

get_nodeheights <- function(tree){
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
    nh <- nh - min(nh)
    tiph <- nh[idx_tip]
    ind <- edge[idx_tip, 2]
    names(tiph) <- tpo$tip.label[ind]
    list(nodeheights=nh[idx_internal], tipheights=tiph, nh=nh)
}
