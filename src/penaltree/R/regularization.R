
get_gpnet <- function(x, y, calc_convex_obj, get_scale, calcalpha=1, nlambda=100,
                      lambda.min.ratio=0.01, lambda=NULL, standardize=TRUE,
                      intercept=TRUE, thresh=1e-4, dfmax=nvars + 1,
                      pmax=min(dfmax*2 + 20,nvars), exclude,
                      penalty.factor=rep(1, nvars), lower.limits=-Inf,
                      upper.limits=Inf, maxit=100){
    if (alpha > 1) {
        warning("alpha >1; set to 1")
        alpha <- 1
    }
    if (alpha < 0) {
        warning("alpha<0; set to 0")
        alpha <- 0
    }
    alpha <- as.double(alpha)
    this.call <- match.call()
    nlam <- as.integer(nlambda)
    np <- dim(x)
    if (is.null(np) | (np[2] < 1))
        stop("x should be a matrix with 1 or more columns")
    nrates <- as.integer(np[1])
    nvars <- as.integer(np[2])
    k <- nchar(y$tmlol[[1]][[1]]$alphabet)
    if (k*(k-1) != nrates)
        stop(paste("number of ordered pairs of states for trait (", k*(k-1),
                   ") not consistent with the number of predicted transition ",
                   "rates (", nrates, ")", sep = ""))
    vnames <- colnames(x)
    if (is.null(vnames))
        vnames <- paste("V", seq(nvars), sep = "")
    ne <- as.integer(dfmax)
    nx <- as.integer(pmax)
    if (!missing(exclude)) {
        jd <- match(exclude, seq(nvars), 0)
        if (!all(jd > 0))
            stop("Some excluded variables out of range")
        jd <- as.integer(c(length(jd), jd))
    }
    else jd <- as.integer(0)
    vp <- as.double(penalty.factor)
    if (any(lower.limits > 0)) {
        stop("Lower limits should be non-positive")
    }
    if (any(upper.limits < 0)) {
        stop("Upper limits should be non-negative")
    }
    if (length(lower.limits) < nvars) {
        if (length(lower.limits) == 1)
            lower.limits <- rep(lower.limits, nvars)
        else stop("Require length 1 or nvars lower.limits")
    }
    else lower.limits <- lower.limits[seq(nvars)]
    if (length(upper.limits) < nvars) {
        if (length(upper.limits) == 1)
            upper.limits = rep(upper.limits, nvars)
        else stop("Require length 1 or nvars upper.limits")
    }
    else upper.limits <- upper.limits[seq(nvars)]
    cl = rbind(lower.limits, upper.limits)
    if (any(abs(cl) < thresh)) {
        stop("Cannot enforce limits this close to zero")
    }
    storage.mode(cl) <- "double"
    isd <- as.integer(standardize)
    intr <- as.integer(intercept)
    thresh <- as.double(thresh)
    if (is.null(lambda)) {
        if (lambda.min.ratio >= 1)
            stop("lambda.min.ratio should be less than 1")
        flmin <- as.double(lambda.min.ratio)
        ulam <- double(1)
    }
    else {
        flmin <- as.double(1)
        if (any(lambda < 0))
            stop("lambdas should be non-negative")
        ulam <- as.double(rev(sort(lambda)))
        nlam <- as.integer(length(lambda))
    }
    fit <- gpnet(x, y, calc_convex_obj, get_scale, alpha, nobs, nvars, jd, vp,
                 cl, ne, nx, nlam, flmin, ulam, thresh, isd, intr, vnames,
                 maxit)
    fit$call <- this.call
    fit$nrates <- nrates
    class(fit) <- c(class(fit), "gpnet")
    fit
}

gpnet <- function(x, y, calc_convex_obj, get_scale, alpha, nobs, nvars, jd, vp,
                  cl, ne, nx, nlam, flmin, ulam, thresh, isd, intr, vnames,
                  maxit, a=0.1, r=0.01, relStart=0.1, mubar=1, beta=0.9,
                  verbose=FALSE, debug=TRUE, initFactor=10){
    maxit <- as.integer(maxit)
    niter <- 0
    if(!intr) stop("not implemented")
    dim <- nvars + intr
    parInds <- seq(dim)
    I <- diag(nrow=dim)
    nll <- function(w){
        -calc_convex_loglik(w=w, x=x, y=y)
    }
    scaleEst <- get_scale(y)
    par <- c(scaleEst, rep(0, nvars))
    gnll <- numDeriv::grad(nll, x=par, method='simple')
    mu <- mubar
    stopifnot(beta>0, beta<1)
    G <- diag(initFactor * abs(gnll), ncol=dim)
    if(flmin<1){
        lstart <- max(abs(gnll))
        loglstart <- log10(lstart) + relStart
        log10LambdaRange <- -log10(flmin)
        loglend <- loglstart - log10LambdaRange
        logstep <- (loglend - loglstart)/nlam
        lambda <- loglstart + 0:(nlam - 1)*logstep
        lambda <- 10^lambda
    } else {
        lambda <- ulam
    }
    res <- list()
    fsg <- function(p, g, l1, l2, h) {
        if(p==0) {
            max(abs(g) - l1, 0)
        } else if (p > 0){
            (-g - l1 - l2*p)/(h + l2)
        } else {
            (-g + l1 - l2*p)/(h + l2)
        }
    }
    for (i in seq_along(lambda)){
        l1penalty <- lambda[i] * c(0, vp) * alpha
        l2penalty <- lambda[i] * c(0, vp) * (1 - alpha)
        penFunc <- function(par) abs(par) %*% l1penalty + (par^2 %*% l2penalty)/2
        k <- 0
        nlp <- nll(par)
        F1 <- nlp + penFunc(par)
        H <- I/(2*mu) + G
        sg <- mapply(fsg, p=par, g=gnll, l1=l1penalty, l2=l2penalty, h=diag(H))
        nsg <- max(abs(sg))
        while (nsg > thresh && k < maxit){
            H <- I/(2*mu) + G
            d <- numeric(dim)
            dlist <- list()
            fmlist <- list()
            inactive <- par != 0 | sg != 0 ## inactive means not actively fixed to zero in solution
            nInactive <- sum(inactive)
            ndesc <- (1 + floor(k * a)) * nInactive
            for (nd in 1:ndesc){
                j <- sample.int(n=nInactive, size=1)
                j <- parInds[inactive][j]
                Hd <- H %*% d
                gr <- gnll[j] + Hd[j]
                if (par[j] + d[j] > 0 || (par[j] + d[j] == 0 & -gr > 0)){
                    z <- (-gr - l1penalty[j] - l2penalty[j]*(par[j] + d[j]))/(H[j,j] + l2penalty[j])
                    if (par[j] + d[j] + z < 0){
                        d[j] <- -par[j]
                    } else {
                        d[j] <- d[j] + z
                    }
                } else {
                    z <- (-gr + l1penalty[j] - l2penalty[j]*(par[j] + d[j]))/(H[j,j] + l2penalty[j])
                    if (par[j] + d[j] + z > 0){
                        d[j] <- -par[j]
                    } else {
                        d[j] <- d[j] + z
                    }
                }
                dlist[[nd]] <- d
                fmlist[[nd]] <- d %*% (H/2) %*% d + gnll %*% d + F1 - penFunc(par) + penFunc(par + d)
            }
            if(any(diff(c(F1, unlist(fmlist)))>1e-10)) browser()
            par2 <- par + d
            if (any(par2[-1] > cl[2, ] || any(par2[-1] < cl[1, ]))){
                if(verbose) cat('backtracking: out of bounds', '\n')
                mu <- mu * beta
            } else {
                Fmod <- d %*% (H/2) %*% d + gnll %*% d + F1 - penFunc(par) + penFunc(par2)
                if(Fmod  > F1 && !isTRUE(all.equal(Fmod, F1))) {
                    if (debug) browser()
                    if(verbose) cat('backtracking: poorly solved model', '\n')
                    mu <- mu * beta
                } else {
                    nlp <- nll(par2)
                    F2 <- nlp + penFunc(par2)
                    if (F2 - F1 > r * (Fmod - F1)){
                        if(verbose) cat('backtracking: insufficient decrease', '\n')
                        mu <- mu * beta
                    } else {
                        gnll2 <- numDeriv::grad(nll, x=par2, method='simple')
                        yvec <- gnll2 - gnll
                        s <- d
                        ys <- yvec %*% s
                        if (ys > 0){
                            ## Hessian approximation update via BFGS via 8.19 in Nocedal and Wright
                            sColVec <- matrix(s, ncol=1)
                            yColVec <- matrix(yvec, ncol=1)
                            M <- G %*% sColVec %*% t(sColVec) %*% G
                            M <- M / as.numeric(t(sColVec) %*% G %*% sColVec)
                            M <- G - M
                            G <- M + yColVec %*% t(yColVec) / as.numeric(t(yColVec) %*% sColVec)
                            if (any(diag(G) <= 0)) browser()
                        } else if(verbose){
                            cat('skipping Hessian update: ys <= 0', '\n')
                        }
                        ## update vars
                        k <- k + 1
                        par <- par2
                        gnll <- gnll2
                        dF <- F2 - F1
                        F1 <- F2
                        H <- I/(2*mu) + G
                        sg <- mapply(fsg, p=par, g=gnll, l1=l1penalty, l2=l2penalty, h=diag(H))
                        nsg <- max(abs(sg))
                        if (debug && mu < 1e-8 && nsg > thresh) browser()
                        if (verbose) {
                            cat('lambda: ', lambda[i], '\n')
                            cat('k: ', k, '\n')
                            cat('F: ', as.numeric(F1), '\n')
                            cat('par: ', signif(par, 3), '\n')
                            cat('grad: ', signif(gnll, 3), '\n')
                            cat('nsg: ', signif(nsg, 3), '\n')
                            cat('\n')
                        }
                    }
                }
            }
        }
        convergence <- ifelse(k == maxit, 'no', 'yes')
        res[[i]] <- list(par=par, nll=nlp, k=k, lambda=lambda[i], convergence=convergence, mu=mu, nsg=nsg, sg=sg)
        mu <- mubar
    }
    path <- sapply(res, '[[', 'par')
    ret <- list(a0=path[1,])
    beta <- path[-1, ]
    colnames(beta) <- paste("s", seq(ncol(beta)) - 1, sep = "")
    rownames(beta) <- vnames
    ret$beta <- beta
    ret$lambda <- sapply(res, '[[', 'lambda')
    ret$nll <- sapply(res, '[[', 'nll')
    ret$df <- colSums(abs(beta) > 0)
    ret$dim <- dim(x)
    ret$niterations <- sum(sapply(res, '[[', 'k'))
    ret$jerr <- paste('convergence: ', sapply(res, '[[', 'convergence'))
    if(verbose) cat('Completed regularization path', '\n')
    ret
}
