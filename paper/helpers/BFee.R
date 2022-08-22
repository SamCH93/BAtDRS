## function to compute equality of effect size BF
BFee_ <- function(to, tr, so, sr, tau, Q = (to - tr)^2/(so^2 + sr^2),
                  H = tau^2/((so^2 + sr^2)/2)) {
    ## check inputs
    stopifnot(
        ## length(to) == 1,
        ## is.numeric(to),
        ## is.finite(to),

        ## length(tr) == 1,
        ## is.numeric(tr),
        ## is.finite(tr),

        ## length(so) == 1,
        ## is.numeric(so),
        ## is.finite(so),
        ## 0 < so,

        ## length(sr) == 1,
        ## is.numeric(sr),
        ## is.finite(sr),
        ## 0 < sr,

        length(Q) == 1,
        is.numeric(Q),
        is.finite(Q),
        0 <= Q,

        length(H) == 1,
        is.numeric(H),
        is.finite(H),
        0 < H
    )

    bf <- sqrt(1 + H) * exp(-0.5 * Q * H / (1 + H))
    return(bf)
}
BFee <- Vectorize(FUN = BFee_)

## function to compute probability of replication success for equality of
## effect size BF
powerBFee_ <- function(level, to, so, sr, tauMT,
                       designPrior = c("predictive", "EB", "conditional", "H0"),
                       tau = 0) {
     ## input checks
    stopifnot(length(level) == 1,
              is.numeric(level),
              is.finite(level),
              0 < level,

              length(to) == 1,
              is.numeric(to),
              is.finite(to),

              length(so) == 1,
              is.numeric(so),
              is.finite(so),
              0 < so,

              length(sr) == 1,
              is.numeric(sr),
              is.finite(sr),
              0 < sr,

              length(tauMT) == 1,
              is.numeric(tauMT),
              is.finite(tauMT),
              0 < tauMT,

              !is.null(designPrior))
    designPrior <- match.arg(designPrior)

    stopifnot(length(tau) == 1,
              is.numeric(tau),
              is.finite(tau),
              0 <= tau)


    ## compute some fixed quantities
    H <- tauMT^2/(so^2 + sr^2)*2
    X2 <- (log(1 + H) - log(level^2))*(1 + 1/H)*(so^2 + sr^2)

    ## check whether possible to achieve power > 0
    if (1 + H < level ^2) {
        pow <- 0
    } else {

        ## determine parameters of predictive distribution of replication estimate
        if (designPrior == "predictive") {
            mu <- to
            sigma2 <- so^2 + sr^2 + 2*tau^2
        } else if (designPrior == "EB") {
            zo <- to/so
            h <- tau^2/so^2
            s <- pmax(1 - (1 + h)/zo^2, 0)
            mu <- s*to
            sigma2 <- s*(so^2 + tau^2) + sr^2 + tau^2
        } else if (designPrior == "H0") {
            mu <- 0
            sigma2 <- sr^2
        } else { ## conditional design prior
            mu <- to
            sigma2 <- sr^2
        }


        ## compute probability of replication success
        pow <- stats::pchisq(q = X2/sigma2, df = 1,
                             ncp = (mu - to)^2/sigma2, lower.tail = TRUE)
    }
    return(pow)
}
powerBFee <- Vectorize(FUN = powerBFee_)

## function to compute the sample size based on the BFee
## TODO write function
ssBFee_ <- function(level, power, to, so, tauMT,
                    designPrior = c("predictive", "EB", "conditional"),
                    tau = 0, searchFac = 1) {
    ## input checks
    stopifnot(length(level) == 1,
              is.numeric(level),
              is.finite(level),
              0 < level,

              length(power) == 1,
              is.numeric(power),
              is.finite(power),
              0 < power, power < 1,

              length(to) == 1,
              is.numeric(to),
              is.finite(to),

              length(so) == 1,
              is.numeric(so),
              is.finite(so),
              0 < so,

              length(tauMT) == 1,
              is.numeric(tauMT),
              is.finite(tauMT),
              0 < tauMT,

              !is.null(designPrior))
    designPrior <- match.arg(designPrior)

    stopifnot(length(tau) == 1,
              is.numeric(tau),
              is.finite(tau),
              0 <= tau,

              length(searchFac) == 1,
              is.numeric(searchFac),
              is.finite(searchFac),
              0 < searchFac)


    ## function to search with uniroot over log relative variance
    rootFun <- function(logc) {
        sr <- so^2/exp(logc)
            suppressWarnings({
                res <- powerBFee(level = level, to = to, so = so, sr = sr,
                                 tauMT = tauMT, designPrior = designPrior,
                                 tau = tau) - power
            })
            return(res)
    }
    res <- try(stats::uniroot(f = rootFun, interval = c(-100, 100)*searchFac)$root,
               silent = TRUE)
        if (class(res) == "try-error") {
            sr <- NaN
            t1e <- NaN
        } else {
            sr <- so^2/exp(res)
        }

     ## return relative variance, actual error rates, and inputs
    if (!is.na(sr)) {
        suppressWarnings({
            power <- powerBFee(level = level, to = to, so = so, sr = sr,
                               tauMT = tauMT, designPrior = designPrior,
                               tau = tau)
            t1e <- powerBFee(level = level, to = to, so = so, sr = sr,
                             tauMT = tauMT, designPrior = "H0",
                             tau = tau)
        })
    }
    out <- list(sr = sr, level = level, t1e = t1e, power = power, to = to,
                so = so, tauMT = tauMT, designPrior = designPrior, tau = tau)
    return(out)
}

ssBFee__ <- Vectorize(FUN = ssBFee_)
ssBFee <- function(level, power, to, so, tauMT, designPrior = "predictive",
                   tau = 0, searchFac = 1) {
    t(ssBFee__(level, power, to, so, tauMT, designPrior = designPrior,
                   tau = tau, searchFac = searchFac))
}
