#' Generate simulated data
#'
#' This function generates simulated data.
#'
#' @param seed Random seed (optional)
#' @param p Variable dimension, default is 10
#' @param ni Sample size per group, default is 100
#' @param K Number of groups, default is 10
#' @param error.Dis Error distribution, default is 'n'
#' @param tau Quantile level, default is 1/2
#' @param dgp Data generating process, default is 'dgp1'
#' @param ... Other parameters
#'
#' @return Returns a list containing simulated data including Y.hat (predicted values),
#'   X (independent variables), tB (regression coefficient matrix), Y.obs (observed values),
#'   and para.input (input parameters).
#'
#' @examples
#' dgp <- 'dgp1'
#' p <- 400
#' ni <- 20
#' K <- 10
#' error_dis <- 't'
#' tau <- 3/4
#' intercept <- !tau == 1/2
#' para.cases <- list(dgp = dgp, p = p, ni = ni, K = K, error.Dis = error_dis, tau = tau) # para for cases
#' data <- datasim(seed = 0, p = p, ni = ni, K = K, error.Dis = error_dis, tau = tau, dgp = dgp)
#'
#' @export
datasim<- function(
        seed = NULL,
        p = 10,
        ni = 100,
        K = 10,
        error.Dis = 'n',
        tau = 1 / 2,
        dgp = 'dgp1',
        ...
) {
    dgp <- as.character(dgp)
    tau <- as.numeric(tau)
    intercept <- !tau == 1/2
    para.input <- list(seed = seed, p = p, ni = ni, K = K, error.Dis = error.Dis, tau = tau, dgp = dgp, intercept = intercept)
    if (!is.null(seed)) {
        set.seed(seed)
    }

    if (dgp == 'dgp1') {
        beta1 <- c(rep(-1, 5), rep(1, 5))
        beta2 <- c(rep(2, 5), rep(2, 5))
        beta <- matrix(data = c(beta1, beta2, rep(0, K * p - 2 * K)), K, p, byrow = FALSE)
    } else if (dgp == 'dgp2') {
        beta1 <- c(rep(-2, 3), rep(0, 4), rep(2, 3))
        beta2 <- rep(1, 10)
        beta3 <- c(rep(-2, 4), rep(2, 4), rep(0, 2))
        beta <- matrix(data = c(beta1, beta2, beta3, rep(0, K * p - 3 * K)), K, p, byrow = FALSE)
    }

    X.temp <- NULL
    Y.hat <- NULL
    tB <- c(NULL)

    if (intercept) {
        if (error.Dis == 't') {
            beta0 <- rep(qt(tau, df = 2), K)
        } else if (error.Dis == 'n') {
            beta0 <- rep(qnorm(tau), K)
        } else if (error.Dis == 'mix') {
            b <- rbinom(1e6, size = 1, 0.9)
            bb <- b * rnorm(1e6) + (1 - b) * rcauchy(1e6)
            beta0 <- rep(quantile(bb, tau), K)
            names(beta0) <- NULL
        }
    }

    for (k in 1:K) {
        Cx <- 0.7 * diag(1, p) + 0.3 * matrix(1, p, p)
        x <- matrix(rnorm(ni * p, 0, 1), nrow = ni, ncol = p)
        X.temp[[k]] <- x %*% chol(Cx)
        Y.hat[[k]] <- c(X.temp[[k]] %*% beta[k,])

        if (intercept) {
            tB <- cbind(tB, c(beta0[k], beta[k,]))
        } else {
            tB <- cbind(tB, beta[k,])
        }
    }

    if (intercept) {
        X.temp <- lapply(X.temp, function(arr) {
            cbind(rep(1, dim(arr)[1]), arr)
        })
    }

    para.input <- list(
        seed = seed,
        p = p,
        ni = ni,
        K = K,
        error.Dis = error.Dis,
        tau = tau,
        dgp = dgp,
        intercept = intercept


    )

    Y.temp <- NULL

    for (k in 1:K) {
        if (!is.null(seed)) {
            set.seed(seed)
        }

        if (error.Dis == 't') {
            e1 <- rt(n = ni, df = 2)
        } else if (error.Dis == 'n') {
            e1 <- rnorm(ni)
        } else {
            b <- rbinom(ni, size = 1, 0.9)
            e1 <- b * rnorm(ni) + (1 - b) * rcauchy(ni)
        }

        e <- as.vector(t(e1))
        Y.temp[[k]] <- Y.hat[[k]] + e
    }

    list(
        Y.hat = Y.hat,
        X = X.temp,
        tB = tB,
        Y.obs = Y.temp,
        para.input = para.input
    )
}
