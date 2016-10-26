
set_branchlengths <- function(tree, nodeheights){
    tree <- reorder(tree, order="postorder")
    edge <- tree$edge
    stopifnot(nrow(edge) == length(nodeheights) - 1)
    ipheight <- match(edge[, 1], edge[,2])
    ipheight[is.na(ipheight)] <- nrow(edge) + 1
    ## Ensure ordering of nodeheights is valid
    for (i in seq(1, nrow(edge))){
        if (nodeheights[ipheight[i]] < nodeheights[i]){
            nodeheights[ipheight[i]] <- nodeheights[i] 
        }
    }
    ## Set branch lengths to differences in heights
    for (i in seq(1, nrow(edge))){
        tree$edge.length[i] <- nodeheights[ipheight[i]] - nodeheights[i]
    }
    tree
}

pml_wrapper <- function(tree, nodeheights, phydata, ...){
    tree <- set_branchlengths(tree, nodeheights)
    pml(tree, data=phydata, ...)$logLik
}
