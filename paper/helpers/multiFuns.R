powerBFsmulti_ <- function(level, zo, c,
                           designPrior = c("predictive", "H0"),
                           h = 0, mu = 0, n = 1) {

    ## input checks
    stopifnot(length(level) == 1,
              is.numeric(level),
              is.finite(level),
              0 < level, level <= 1,

              length(zo) == 1,
              is.numeric(zo),
              is.finite(zo),

              length(c) == 1,
              is.numeric(c),
              is.finite(c),
              0 < c,

              !is.null(designPrior))
    designPrior <- match.arg(designPrior)

    stopifnot(length(h) == 1,
              is.numeric(h),
              is.finite(h),
              0 <= h,

              length(mu) == 1,
              is.numeric(mu),
              is.finite(mu),

              length(n) == 1,
              is.numeric(n),
              is.finite(n),
              0 < n)


    ## determine parameters of predictive distribution of zr
    if (designPrior == "predictive") {
        mu <- zo*sqrt(c)
        sigma2 <- 1 + n*(1 + h)/(h + n/c)
    } else if (designPrior == "H0") {
        mu <- 0
        sigma2 <- 1
    }
    ## else if (designPrior == "EB") {
    ##     s <- pmax(1 - (1 + h)/zo^2, 0)
    ##     mu <- s*zo*sqrt(n*c)
    ##     sigma2 <- (1 + c*h)/n + s*c*(1 + h)
    ## } else { ## conditional design prior
    ##     ## TODO implement for multisite
    ##     mu <- zo*sqrt(c)
    ##     sigma2 <- 1
    ## }

    ## compute sufficiently sceptical prior variance for specified level
    g <- vss(x = zo, gamma = level)
    if (is.nan(g)) {
        pow <- 0
    } else {
        ## compute probability that RS at level
        ## TODO: find solution for numerical problems with non-central chi-squared
        ## calculations when g is a bit above or below 1, currently set numerical
        ## tolerance quite large
        if (isTRUE(all.equal(g, 1, tolerance = 1e-2))) {
            D <- (-log(2) + zo^2*(0.5 + 1/(1/c + 1)))/(2*zo*sqrt(c))*(1 + c)
            ## D <- 0.5*(zo*sqrt(c) + (1/c + 1)*sqrt(c)*(0.5*zo - log(2)/zo))
            pow <- stats::pnorm(q = sign(zo)*(mu - D)/sqrt(sigma2), mean = 0, sd = 1)
            ## if (zo > 0) pow <- stats::pnorm(q = D, mean = mu, sd = sqrt(sigma2),
            ##                                 lower.tail = FALSE)
            ## else pow <- stats::pnorm(q = D, mean = mu, sd = sqrt(sigma2))
        } else {
            A <- log((1/c + 1)/(1/c + g)/(1 + g)) + zo^2/(1/g + 1) + zo^2/(1 -g)
            B <- c*(1 - g)/(1 + c*g)/(1 + c)
            M <- zo*(1 + c*g)/sqrt(c)/(g - 1)
            lambda <- (mu - M)^2/sigma2
            res <- stats::pchisq(q = A/B/sigma2, df = 1, ncp = lambda,
                                 lower.tail = FALSE)
            ## when g > 1, the inequality flips and we need to take the complement...
            if (g < 1) pow <- res
            else pow <- 1 - res
        }
    }
    return(pow)
}
powerBFsmulti <- Vectorize(FUN = powerBFsmulti_)

ssBFsmulti_ <- function(level = NA, t1e = NA, power, zo, h = 0,
                        designPrior = c("predictive"),
                        n = 1,
                        searchFac = 1) {

    ## input checks
    stopifnot(length(level) == 1,
              length(t1e) == 1,
              xor(is.na(t1e), is.na(level)),
              xor(is.numeric(t1e), is.numeric(level)),

              length(power) == 1,
              is.numeric(power),
              is.finite(power),
              0 < power, power < 1,

              length(zo) == 1,
              is.numeric(zo),
              is.finite(zo),

              length(h) == 1,
              is.numeric(h),
              is.finite(h),
              0 <= h,

              !is.null(designPrior))
    designPrior <- match.arg(designPrior)

    stopifnot(length(n) == 1,
              is.numeric(n),
              is.finite(n),
              0 < n,

              length(searchFac) == 1,
              is.finite(searchFac),
              0 < searchFac)

    ## absolut cut-off approach
    if (!is.na(level)) {
        stopifnot(is.finite(level),
                  0 < level, level < 1)

        ## function to search with uniroot over relative variance c
        rootFun <- function(logc) {
            res <- powerBFsmulti(level = level, zo = zo, c = exp(logc), h = h,
                                 designPrior = designPrior, n = n) - power
            return(res)
        }
        res <- try(stats::uniroot(f = rootFun, interval = c(-5, 5)*searchFac)$root,
                   silent = TRUE)
        if (class(res) == "try-error") {
            c <- NaN
        } else {
            c <- exp(res)
        }

    ## relative cut-off approach
    } else {
        stopifnot(is.finite(t1e),
                  0 < t1e, t1e < power)

        ## function to search with uniroot over level for fixed t1e and c
        levelFun <- function(loglevel, c) {
            res <- powerBFsmulti(level = exp(loglevel), zo = zo, c = c, h = h,
                                 designPrior = "H0", n = n) - t1e
            return(res)
        }
        ## function to search with uniroot over relative variance c
        rootFun <- function(logc) {
            loglevel <- try(stats::uniroot(f = levelFun,
                                           interval = c(-5, 0)*searchFac,
                                           c = exp(logc))$root, silent = TRUE)
            if (class(loglevel) == "try-error") res <- NaN
            else {
                res <- powerBFsmulti(level = exp(loglevel), zo = zo,
                                     c = exp(logc), h = h,
                                     designPrior = designPrior, n = n) - power
            }
            return(res)
        }
        ## HACK evaluate root function for some logc to find limits
        ## for numerical search which are not NaN
        logcTest <- seq(-3, 3, length.out = 20)*searchFac
        rootTest <- sapply(X = logcTest, rootFun)
        low <- logcTest[utils::tail(which(!is.na(rootTest) & rootTest < 0), n = 1)]
        up <- logcTest[utils::head(which(!is.na(rootTest) & rootTest > 0), n = 1)]

        res <- try(stats::uniroot(f = rootFun, interval = c(low, up))$root, silent = TRUE)
        if (class(res) == "try-error") {
            c <- NaN
            t1e <- NaN
        } else {
            c <- exp(res)
            level <- exp(stats::uniroot(f = levelFun, interval = c(-5, 0)*searchFac,
                                        c = c)$root)
        }
    }

    ## return relative variance, actual error rates, and inputs
    if (!is.na(c)) {
        power <- powerBFsmulti(level = level, zo = zo, c = c, h = h,
                               designPrior = designPrior, n = n)
        t1e <- powerBFsmulti(level = level, zo = zo, c = c, h = h,
                             designPrior = "H0", n = n)
    }
    out <- list(c = c, level = level, t1e = t1e, power = power, zo = zo, h = h,
                designPrior = designPrior, n = n)
    return(out)
}

## stack in a matrix to transpose afterwards
ssBFsmulti__ <- Vectorize(FUN = ssBFsmulti_, SIMPLIFY = TRUE)
ssBFsmulti <- function(level = NA, t1e = NA, power, zo, h = 0,
                  designPrior = "predictive", n = 1, searchFac = 1) {
    t(ssBFsmulti__(level = level, t1e = t1e, power = power, zo = zo, h = h,
                   designPrior = designPrior, n = n, searchFac = searchFac))
}


powerBFrmulti_ <- function(level, zo, c, g = 0,
                           designPrior = c("predictive", "H0"),
                           h = 0, mu = 0, n = 1) {
    ## input checks
    stopifnot(length(level) == 1,
              is.numeric(level),
              is.finite(level),
              0 < level,

              length(zo) == 1,
              is.numeric(zo),
              is.finite(zo),

              length(c) == 1,
              is.numeric(c),
              is.finite(c),
              0 <= c,

              length(g) == 1,
              is.numeric(g),
              is.finite(g),
              0 <= g,

              !is.null(designPrior))
    designPrior <- match.arg(designPrior)

    stopifnot(length(h) == 1,
              is.numeric(h),
              is.finite(h),
              0 <= h,

              length(mu) == 1,
              is.numeric(mu),
              is.finite(mu),

              length(n) == 1,
              is.numeric(n),
              is.finite(n),
              0 < n)

    ## determine parameters of predictive distribution of zr
    if (designPrior == "predictive") {
        mu <- zo*sqrt(c)
        sigma2 <- 1 + n*(1 + h)/(h + n/c)
    } else if (designPrior == "H0") {
        mu <- 0
        sigma2 <- 1
    }
    ## else if (designPrior == "EB") {
    ##     s <- pmax(1 - (1 + h)/zo^2, 0)
    ##     mu <- s*zo*sqrt(n*c)
    ##     sigma2 <- (1 + c*h)/n + s*c*(1 + h)
    ## } else { ## conditional design prior
    ##     ## TODO implement for multisite
    ##     mu <- zo*sqrt(c)
    ##     sigma2 <- 1
    ## }


    ## compute probability that BFr < level
    k <- log(level)
    if (isTRUE(all.equal(g, 1))) {
        D <- (k + 0.5*zo^2/(1/c + 1))*(1 + c)/(zo*sqrt(c))
        pow <- stats::pnorm(q = sign(zo)*(D - mu)/sqrt(sigma2), mean = 0, sd = 1)
    } else {
        A <- log((1 + c)/(1 + g*c)) - 2*k + zo^2/(1 - g)
        B <- (1 - g)/(1 + c*g)/(1 + 1/c)
        M <- zo*sqrt(c)*(1/c + g)/(g - 1)
        lambda <- (mu - M)^2/sigma2
        prob <- stats::pchisq(q = A/B/sigma2, df = 1, ncp = lambda,
                              lower.tail = FALSE)

        ## when g > 1, the inequality flips and we need to take the complement
        if (g > 1) pow <- 1 - prob
        else pow <- prob
    }
    return(pow)
}
powerBFrmulti <- Vectorize(FUN = powerBFrmulti_)

ssBFrmulti_ <- function(level = NA, t1e = NA, power, zo, h = 0, g = 0,
                        designPrior = c("predictive"),
                        n = 1, searchFac = 1) {
    ## input checks
    stopifnot(length(level) == 1,
              length(t1e) == 1,
              xor(is.na(t1e), is.na(level)),
              xor(is.numeric(t1e), is.numeric(level)),

              length(power) == 1,
              is.numeric(power),
              is.finite(power),
              0 < power, power < 1,

              length(zo) == 1,
              is.numeric(zo),
              is.finite(zo),

              length(h) == 1,
              is.numeric(h),
              is.finite(h),
              0 <= h,

              length(g) == 1,
              is.numeric(g),
              is.finite(g),
              0 <= g,

              !is.null(designPrior))
    designPrior <- match.arg(designPrior)

    stopifnot(length(n) == 1,
              is.numeric(n),
              is.finite(n),
              0 < n,

              length(searchFac) == 1,
              is.finite(searchFac),
              0 < searchFac)

    ## absolut cut-off approach
    if (!is.na(level)) {
        stopifnot(is.finite(level),
                  0 < level, level < 1)

        ## function to search with uniroot over relative variance c
        rootFun <- function(logc) {
            suppressWarnings({
                res <- powerBFrmulti(level = level, zo = zo, c = exp(logc),
                                     h = h, g = g, designPrior = designPrior,
                                     n = n) - power
            })
            return(res)
        }
        res <- try(stats::uniroot(f = rootFun, interval = c(-5, 5)*searchFac)$root,
                   silent = TRUE)
        if (class(res) == "try-error") {
            c <- NaN
            t1e <- NaN
        } else {
            c <- exp(res)
        }

    ## relative cut-off approach
    } else {
        stopifnot(is.finite(t1e),
                  0 < t1e, t1e < power)

        ## function to search with uniroot over level for fixed t1e and c
        levelFun <- function(loglevel, c) {
            suppressWarnings({
                res <- powerBFrmulti(level = exp(loglevel), zo = zo, c = c,
                                     h = h, g = g, designPrior = "H0", n = n) - t1e
            })
            return(res)
        }
        ## function to search with uniroot over relative variance c
        rootFun <- function(logc) {
            suppressWarnings({
                loglevel <- try(stats::uniroot(f = levelFun,
                                               interval = c(-5, 0)*searchFac,
                                               c = exp(logc))$root, silent = TRUE)
                if (class(loglevel) == "try-error") res <- NaN
                else {
                    res <- powerBFrmulti(level = exp(loglevel), zo = zo,
                                         c = exp(logc), h = h, g = g,
                                         designPrior = designPrior, n = n) - power
                }
            })
            return(res)
        }
        ## HACK evaluate root function for some logc to find limits
        ## for numerical search which are not NaN
        logcTest <- seq(-3, 3, length.out = 20)*searchFac
        rootTest <- sapply(X = logcTest, rootFun)
        low <- logcTest[utils::tail(which(!is.na(rootTest) & rootTest < 0), n = 1)]
        up <- logcTest[utils::head(which(!is.na(rootTest) & rootTest > 0), n = 1)]
        res <- try(stats::uniroot(f = rootFun, interval = c(low, up))$root,
                   silent = TRUE)
        if (class(res) == "try-error") {
            c <- NaN
        } else {
            c <- exp(res)
            level <- exp(stats::uniroot(f = levelFun, interval = c(-5, 0)*searchFac,
                                        c = c)$root)
        }
    }

    ## return relative variance, actual error rates, and inputs
    if (!is.na(c)) {
        suppressWarnings({
            power <- powerBFrmulti(level = level, zo = zo, c = c, h = h, g = g,
                              designPrior = designPrior, n = n)
            t1e <- powerBFrmulti(level = level, zo = zo, c = c, h = h, g = g,
                                 designPrior = "H0", n = n)
        })
    }
    out <- list(c = c, level = level, t1e = t1e, power = power, zo = zo, h = h,
                g = g, designPrior = designPrior, n = n)
    return(out)
}

## stack in a matrix to transpose afterwards
ssBFrmulti__ <- Vectorize(FUN = ssBFrmulti_, SIMPLIFY = TRUE)
ssBFrmulti <- function(level = NA, t1e = NA, power, zo, h = 0, g = 0,
                       designPrior = "predictive", n = 1, searchFac = 1) {
    t(ssBFrmulti__(level = level, t1e = t1e, power = power, zo = zo, g = g, h = h,
                   designPrior = designPrior, n = n, searchFac = searchFac))
}

## ssBFsmulti(level = 1/3, power = 0.8, zo = 4, n = 5)
## ssBFrmulti(level = 1/3, power = 0.8, zo = 4, n = 5)

## library(BayesRep)
## plot(powerBFsmulti(level = 1/3, zo = 2.5, c = 2, n = seq(1,
## 100), h = 0, designPrior = "predictive"), type = "l")
