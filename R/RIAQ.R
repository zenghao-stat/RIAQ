# TOOL
bic <- function(Y_temp, X_temp, Bmat, type = 'zhang',tau = 1/2,cc=1/10) { #
    # Bmat p x K
    K = length(Y_temp)
    p = dim(X_temp[[1]])[2]
    n = dim(X_temp[[1]])[1]

    # AMS = function(B){
    #     ams = NULL
    #     for (j in 1:ncol(B)) {
    #         ams = c(ams, sum(unique(B[,j]) !=0 ))
    #     }
    #     ams = sum(ams)
    #     return(ams)
    # }

    AMS <-  function(B){ #
        ams = NULL
        B = matrix(B)
        ams = sum(unique(B)!=0 )
        return(ams)
    }

    if (type == 'yang') {
        loss = 0
        for (k in 1:K) {
            loss = loss + norm(Y_temp[[k]] - X_temp[[k]] %*% Bmat[,k],"2")^2
        }
        qsum = AMS((Bmat))
        logmse = log(loss/n)
        logsparsity = cc*log(log(K*p))*(log(n)/n)*qsum #yang 2019
        bic = logmse + logsparsity #Yang 2019
    }

    if (type == 'zhang') {

        rho_tau  <-  function(x,tau){x*(tau - as.numeric(x<0))} #

        loss = 0
        # cc = 8
        for (k in 1:K) {
            loss = loss + sum(rho_tau(Y_temp[[k]] - X_temp[[k]] %*% Bmat[,k],tau))
        }
        qsum = AMS((Bmat))
        logmse = log(loss/n)
        # logsparsity = cc*log(log(n))*(log(K*p)/n)*qsum #Zhang 2019
        logsparsity = cc*log(log(K*p))*(log(n)/n)*qsum #yang 2019
        bic = logmse + logsparsity
    }
    return(list(bic = bic,logmse = logmse, logsparsity= logsparsity ))
}

## Derivative
scad_deriv <- function(x, lambda=1, a=3.7){
    absx <- u <- abs(x)
    u[] <- 0
    index <- absx < a*lambda & absx > 0
    u[ index ] <-
        ifelse( absx[ index ] <= lambda,
                lambda,
                ( a*lambda - absx[ index ] )/( a-1 )
        )
    u[index] <- u[index]*sign( x[index] )
    u[ x == 0 ] <- lambda # because we take derivative as x approaches 0 from above

    u
}


mcp_deriv <- function(x, lambda=1, a=3){
    u <- x
    u[] <- 0
    index <- abs(x) < a*lambda
    u[ index ] <- ifelse( x[index] == 0, lambda, lambda*sign(x[index]) - x[index]/a)
    u
}

# Thresh
thresh_est <- function(z, lambda, tau, a = 3, penalty = c("MCP", "SCAD", "lasso")) {
    ### function to implement the soft-, MCP, SCAD thresholding rule Input
    ### variables:
    ### z: argument
    ### type: thresholding rule 1 = (Adaptive) LASSO
    ### (default) 2 = MCP 3 = SCAD lambda: thresholding level a: default
    ### choice for SCAD penalty.
    ### argmin_x (tau/2)*(x-z)^2 + p_a (|x|,lambda)

    penalty <- match.arg(penalty)

    if (penalty == "lasso") {
        return(
            sign(z) * (abs(z) >= lambda/tau) * (abs(z) - lambda/tau)
        )
    } #soft-threshold(z,lambda/tau)

    if (penalty == "MCP") {
        # a = 2.7 + 1/tau
        return(
            sign(z) *(abs(z) >= lambda/tau) * (abs(z) - lambda/tau)/(1 - 1/(a * tau)) *
                (abs(z) <= a * lambda) +
                (abs(z) > a * lambda)* z
        )
    }

    if (penalty == "SCAD") {
        # a = 3.7 + 1/tau
        return(sign(z) * (abs(z) >= lambda/tau) * (abs(z) - lambda/tau) *(abs(z) <= lambda + lambda/tau) +
                   sign(z) * (abs(z) >= a * lambda/(tau * (a - 1))) *
                   (abs(z) - a * lambda/(tau * (a - 1)))/(1 - 1/((a - 1) * tau)) *
                   (lambda + lambda/tau < abs(z)) *(abs(z) <= a * lambda) +
                   z * (abs(z) > a * lambda))
    }
}


#' RIAQ: Robust Integrative Analysis via Quantile Regression with Homogeneity and Sparsity
#'
#' Implementation of the method proposed in the paper "Robust Integrative Analysis via Quantile Regression with Homogeneity and Sparsity."
#'
#' @param X_temp The input data matrix.
#' @param Y_temp The response variable.
#' @param beta0 Initial value for beta (optional).
#' @param tau The quantile level.
#' @param intercept Logical value indicating whether to include an intercept term.
#' @param pen The penalty function type (e.g., 'MCP', 'SCAD').
#' @param cc A tuning parameter for the penalty function.
#' @param ini.lambda The initial value for the tuning parameter lambda.
#' @param rho1 A tuning parameter.
#' @param rho2 A tuning parameter.
#' @param rho3 A tuning parameter.
#' @param tol Tolerance level for convergence.
#' @param nADMM Maximum number of ADMM iterations.
#' @param is.show Logical value indicating whether to display progress messages.
#' @param para.turning List of turning parameters (optional).
#' @param ... Additional arguments.
#'
#' @return A list containing the estimated coefficients (Bhat).
#'
#' @references
#' Zeng, Hao, Chuang Wan, Wei Zhong, and Tuo Liu. 2023. “Robust Integrative Analysis via Quantile Regression with Homogeneity and Sparsity.” Working paper.
#'
#' @examples
#' \dontrun{
#' dgp = 'dgp1'
#' p = 400
#' ni = 20
#' K = 10
#' error_dis = 't'
#' tau = 3/4
#' intercept = !tau ==1/2
#'
#' para.cases = list(dgp = dgp, p = p, ni = ni, K = K, error.Dis = error_dis, tau = tau) # para for cases
#' data = datasim(seed = 0, p = p, ni = ni, K = K, error.Dis = error_dis,tau = tau,dgp = dgp)
#'
#' para.algo = list(rho1 = 1, rho2 = 1, tol = 1e-4, intercept = intercept , nADMM = 500, pen = 'MCP', cc = 1, ini.lambda = 0.005, is.show = T,para.turning = NULL) # paras for algo
#' args = c(para.cases, para.algo,list(X_temp = data$X, Y_temp = data$Y.obs, beta0 = t(data$tB)))
#' fit= do.call(RIAQ, args = args)
#'
#' Beta_hat = fit$Bhat
#' }
#'
#' @export
RIAQ = function(X_temp, Y_temp, beta0 = NULL, tau, intercept = F, pen = 'MCP', cc = 1, ini.lambda = 0.01, rho1 = 1, rho2 = 1, rho3 = 1, tol = 1e-4, nADMM = 500, is.show = T, para.turning = NULL, ... ) {

    para.input = list(tau = tau, intercept = intercept, pen = pen, cc = round(cc,3), ini.lambda = ini.lambda, rho1 = rho1, rho2 = rho2, rho3 = rho3, tol = tol, nADMM = nADMM, is.show = is.show, para.turning = para.turning )

    cat('paras for algo: ')
    cat(paste(names(para.input),'=', para.input), "\n")

    require(quantreg,quietly = T)

    # Data and common matrix
    {
        p = dim(X_temp[[1]])[2]
        K <- length(X_temp)
        nk <- do.call("cbind",sapply(1:K,function(j) length(Y_temp[[j]]),simplify=F))
        nk <- drop(nk)
        n <- sum(nk)
        E <- NULL
        for (j in 1:(K - 1)) {
            E[[j]] <- matrix(0, nrow = K - j, ncol = K)
            E[[j]][, j] <- 1
            for (k in (j + 1):K) {
                E[[j]][k - j, k] <- -1
            }
        }
        E <- do.call("rbind", E)
        X <- do.call("rbind",X_temp)
        Y <- as.vector(do.call("c",Y_temp))
        XL <- Matrix::bdiag(X_temp)
        D <- (diag(K)+rho1*matrix(1,K,K)/rho2)/(K*rho1+rho2)
        DL <- NULL
        for(i in 1:K){
            Q0 <- sapply(1:K,function(j,Drow,n) Drow[j]*rep(1,n[i])%*%t(rep(1,n[j])), Drow=D[i,],n=nk,simplify = F)
            DL[[i]] <- do.call("cbind",Q0)
        }
        DL <- do.call("rbind",DL)
        AL <- solve(diag(n)+rho3*matrixcalc::hadamard.prod(X%*%t(X),DL))
    }

    if (is.null(para.turning)){
        lambda = 0
        gamma = 0
    }else{
        lambda = para.turning[[1]]
        gamma = para.turning[[2]]
    }
    estimationonce = function(K, nk, beta0 = NULL, X_temp, Y_temp, pen = 'SCAD', ini.lambda, nADMM, Y, tau, intercept = F, rho3, X, rho1, E, rho2, AL, D, gamma, lambda, n, tol , p, show)
    {
        require(rqPen)
        I_temp <- sapply(1:K, function(j, nk) t(rep(1, nk[j])), nk = nk, simplify = F)
        IL <- as.matrix(Matrix::bdiag(I_temp))
        if (is.null(beta0)) {
            beta_temp <- NULL
            z_temp <- NULL
            for (k in 1:K) {
                if (intercept) {
                    ini.rq = rq.lasso.fit(
                        X_temp[[k]][, 2:p],
                        Y_temp[[k]],
                        tau = tau,
                        lambda = ini.lambda,
                        weights = NULL,
                        intercept = T,
                        coef.cutoff = .00001,
                        method = "br"
                    )
                } else{
                    ini.rq <-
                        rq.lasso.fit(
                            X_temp[[k]],
                            Y_temp[[k]],
                            tau = tau,
                            lambda = ini.lambda,
                            weights = NULL,
                            intercept = FALSE,
                            coef.cutoff = .00001,
                            method = "br"
                        )
                }

                beta_temp[[k]] <- round(ini.rq$coef, 3)
                z_temp[[k]] <- ini.rq$residual
                betahat <- do.call("rbind", beta_temp)
                z <- as.vector(unlist(z_temp))
            }
        } else{
            z_temp <- NULL
            for (k in 1:K) {
                z_temp[[k]] = Y_temp[[k]] - X_temp[[k]] %*% beta0[k, ]
            }
            z <- as.vector(unlist(z_temp))
            betahat = beta0
        }

        # cat(betahat)
        alpha_temp <- NULL
        t <- 1
        for (i in 1:(K - 1)) {
            for (j in (i + 1):K) {
                alpha_temp[[t]] <- betahat[i, ] - betahat[j, ]
                t <- t + 1
            }
        }
        alphahat <- do.call("rbind", alpha_temp)
        etahat <- betahat

        anu <- alphahat - alphahat
        bv <- betahat - betahat
        cz <- z - z

        for (iADMM in 1:nADMM) {
            wideY <- Y - z + cz / rho3

            Xy_temp <- sapply(1:length(Y), function(j, mat, list) kronecker(mat[j, , drop = FALSE], list[j]), mat = X, list = wideY, simplify = F)

            Xy = do.call("rbind", Xy_temp)

            Bm <-rho3 * IL %*% Xy + rho1 * t(E) %*% (alphahat - anu / rho1) + rho2 * etahat + bv

            #BLR <- as.vector(AL%*% apply(t(IL)%*%D%*%Bm,1,sum))
            #XIDB <- apply(t(IL))

            BLR <- as.vector((AL %*% apply(X * ( t(IL) %*% D %*% Bm ), 1, sum)))
            BL_temp <- NULL
            BL_temp <- sapply(1:length(BLR), function(j, mat, list)
                kronecker(mat[j, , drop = FALSE], list[j]),
                mat = X, list = BLR, simplify = F)
            BL = do.call("rbind", BL_temp) * rho3  #X is a design matrix after transformation

            betahat <- D %*% (Bm - IL %*% BL)
            # betahat = round(betahat, 4)
            zeta <- betahat - bv / rho2
            xi_temp <- NULL
            t <- 1
            for (i in 1:(K - 1)) {
                for (j in (i + 1):K) {
                    xi_temp[[t]] <- betahat[i, ] - betahat[j, ]
                    t <- t + 1
                }
            }
            xi <- do.call("rbind", xi_temp) + anu / rho1

            etahat <- thresh_est(zeta, gamma, tau = rho2, a = 3, penalty = pen)
            alphahat <- thresh_est(xi, lambda, tau = rho1, a = 3, penalty = pen)

            ##########update z

            z.old <- z

            h_temp <- NULL
            h_temp <- lapply(1:K, function(j) Y_temp[[j]] - X_temp[[j]] %*% betahat[j, ])
            h <- do.call('c', h_temp)

            tmp1 <- h + cz / rho3
            tmp2 <- n * rho3
            z <- rep(0, n)
            loc = tmp1 >= (tau / tmp2)
            z[loc] = (tmp1 - (tau / tmp2))[loc]
            loc = tmp1 <= ((tau - 1) / tmp2)
            z[loc] = (tmp1 - (tau - 1) / tmp2)[loc]

            ####update auxiliary variable

            anu.old <- anu
            bv.old <- bv
            cz.old <- cz

            u1 <- E %*% betahat - alphahat
            u2 <- etahat - betahat
            u3 <- h - z
            anu <- anu.old + rho1 * (u1)
            bv <- bv.old + rho2 * (u2)
            cz <- cz.old + rho3 * (u3)


            primalRes <- norm(u1, type = "F") + norm(u2, type = "F") + norm(u3, type = "2")
            # cat(primalRes)
            if (abs(primalRes) <= tol * p ^ (1 / 2)) {
                if (show) {
                    cat("Convergence !", "\n")
                }

                break
            }

        }
        return(list(Bmat = betahat, iter = iADMM, primalRes = primalRes))
    }
    # ini fit
    ini.fit = estimationonce(K,nk,beta0 = beta0,X_temp, Y_temp, pen = pen, ini.lambda = ini.lambda, nADMM, Y, tau, intercept = intercept, rho3, X,rho1,E,rho2,AL,D,gamma, lambda,n,tol, p,show = F)
    # record ini fit
    {
        beta0 = Bmat = round(ini.fit$Bmat,2)
        iter = ini.fit$iter
        Bmat.uniqueN = apply(Bmat,1, function(ll){length(unique(ll))})
        if (is.show) {cat("0 ini gam = ", round(gamma, 4), "lambda =", round(lambda, 4), 'bic = ', bic(Y_temp, X_temp, t(Bmat), tau = tau, cc = cc)$bic, "\t",Bmat.uniqueN,' iter = ',ini.fit$iter, ' primalRes = ', ini.fit$primalRes, '\n')}
        Blist = list(Bmat) # record all beta.estimated
        path.gamma = c(gamma) # record gamma
        path.lambda = c(lambda) # record lambda
        iterList = c(iter)
        primalResList = c(ini.fit$primalRes)
    }

    # Selection procedure
    if (is.null(para.turning)) {
        if (is.show) {cat("---- begin selection procedure ---- \n")}

        # set turning parameter search path
        {
            paraLen_gamma = 6
            paraLen_lam = 7
            lexp <- - 2
            lexp.delta <- 0.4
            lexp2 <- -3
            lexp.delta2 <- 0.4
            gammalist = round(10 ^ (lexp + seq(1:paraLen_gamma) * lexp.delta), 3)
            lambdalist = round(10 ^ (lexp2 + seq(1:paraLen_lam) * lexp.delta2), 3)
        }



        i = 1
        j = 1
        # while ((!all(Bmat.uniqueN == 1)) & (!all(Bmat.uniqueN == 0))) {
        #     gamma = gammalist[indM[i, 1]]
        #     lambda = lambdalist[indM[i, 2]]
        #     fit = estimationonce( K, nk, beta0 = beta0, X_temp, Y_temp, pen = pen, ini.lambda = ini.lambda, nADMM, Y, tau = tau, intercept = intercept, rho3, X, rho1, E, rho2, AL, D, gamma, lambda, n, tol , p, show = is.show )
        #     Bmat = (round(fit$Bmat, 3))
        #     iter = fit$iter
        #     Bmat.uniqueN = apply(Bmat, 1, function(ll) {length(unique(ll[ll != 0]))})
        #
        #     if (is.show) {
        #         cat( i,' ', "gam = ", round(gamma, 4), "lambda =", round(lambda, 4),'bic = ', bic(Y_temp, X_temp, t(Bmat), tau = tau, cc = cc)$bic, "\t", Bmat.uniqueN, ' iter = ', fit$iter, ' primalRes = ', fit$primalRes,  "\n" )
        #     }
        #
        #     Blist = c(Blist, list(Bmat))
        #     path.gamma = c(path.gamma, gamma)
        #     path.lambda = c(path.lambda, lambda)
        #     iterList = c(iterList, iter)
        #     if (i < paraLen * paraLen) {
        #         i = i + 1
        #     } else{
        #         warning('need larger range of turning parameters')
        #         break
        #     }
        # }


        for (i in 1:length(gammalist)) {
            for (j in 1:length(lambdalist)) {
                gamma = gammalist[i]
                lambda = lambdalist[j]
                fit = estimationonce( K, nk, beta0 = beta0, X_temp, Y_temp, pen = pen, ini.lambda = ini.lambda, nADMM, Y, tau = tau, intercept = intercept, rho3, X, rho1, E, rho2, AL, D, gamma, lambda, n, tol , p, show = is.show )
                Bmat = (round(fit$Bmat, 2))
                iter = fit$iter
                Bmat.uniqueN = apply(Bmat, 1, function(ll) {length(unique(ll[ll != 0]))})

                if (is.show) {
                    cat( i,' ',j, ' ',  "gam = ", round(gamma, 4), "lambda =", round(lambda, 4),'bic = ', bic(Y_temp, X_temp, t(Bmat), tau = tau, cc = cc)$bic, "\t", Bmat.uniqueN, ' iter = ', fit$iter, ' primalRes = ', fit$primalRes,  "\n" )
                }

                Blist = c(Blist, list(Bmat))
                path.gamma = c(path.gamma, gamma)
                path.lambda = c(path.lambda, lambda)
                iterList = c(iterList, iter)
            }
        }

        # bic
        biclist = sapply(Blist, function(ll) {
            bic(Y_temp, X_temp, t(ll), tau = tau, cc = cc)
        })
        return(
            list(
                Bhat = Blist[[which.min(biclist[1, ])]],
                gamma = path.gamma[which.min(biclist[1, ])],
                lambda = path.lambda[which.min(biclist[1, ])],
                iter  = iterList[which.min(biclist[1, ])],
                bic = biclist[1,which.min(biclist[1, ]) ],
                Blist = Blist,
                biclist = biclist,
                path.gamma = path.gamma,
                path.lambda = path.lambda
            )
        )
    } else{
        return(list(Bhat = Bmat, iter = iter,para.input = para.input))
    }
}
