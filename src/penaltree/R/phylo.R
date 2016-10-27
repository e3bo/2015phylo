
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
    nh[idx_tip] <- tipheights
    nh[idx_internal] <- nodeheights
    ipheight <- match(edge[, 1], edge[,2])
    ipheight[is.na(ipheight)] <- nrow(edge) + 1
    ## Ensure ordering of nodeheights is valid
    penalty <-0
    for (i in seq(1, nedge)){
        blen <- nh[ipheight[i]] - nh[i]
        if (blen <= 0){
            penalty <- penalty - blen
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
