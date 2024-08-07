





#' Covariate-Adjusted-Sparse-Matrix-completion
#' Fit function
#'
#' @param y A sparse matrix of class Incomplete.
#' @param X covariate matrix
#' @param svdH (optional) A list consisting of the SVD of the hat matrix. see reduced_hat_decomp function.
#' @param Xterms (optional) A list of terms computed using GetXterms function.
#' @param J (hyperparameter). The maximum rank of the low-rank matrix M = AB. Default is 2
#' @param lambda (hyperparameter) The L2 regularization parameter for A and B
#' @param r (hyperparameter, optional) The rank of covariate effects (beta). The rank will be the same
#'                                      as the covariate matrix if r not provided
#' @return A list of u,d,v of M, Beta, and a vector of the observed Xbeta
#' @examples
#'  CASMC_fit(y,X,J=5)
#' @export
#'
CASMC2_fit <-
  function(y,
           X,
           J = 2,
           r = 2,
           lambda.M = 0,
           lambda.beta = 0,
           # similarity matrix for A
           S.a = NULL,
           lambda.a = 0,
           # similarity matrix for B
           S.b = NULL,
           lambda.b = 0,
           normalized_laplacian = TRUE,
           # convergence parameters
           maxit = 300,
           thresh = 1e-05,
           trace.it = FALSE,
           warm.start = NULL,
           # the following should not be modified
           final.svd = TRUE,
           Qtype = 1,
           learning_rate = 0,
           qiter.max = 10,
           init = "naive",
           min_eigv = 0) {
    stopifnot(inherits(y, "dgCMatrix"))
    irow = y@i
    pcol = y@p
    n <- dim(y)
    m <- n[2]
    n <- n[1]
    k <- ncol(X)
    if (trace.it)
      nz = nnzero(y, na.counted = TRUE)
    #-------------------------------
    laplace.a = laplace.b = F
    if (!is.null(S.a) && lambda.a > 0) {
      laplace.a = T
      L.a = computeLaplacian(S.a, normalized = normalized_laplacian) * lambda.a
    }
    if (!is.null(S.b) && lambda.b > 0) {
      laplace.b = T
      L.b = computeLaplacian(S.b, normalized = normalized_laplacian) * lambda.b
    }
    #---------------------------------------------------
    # warm start or initialize (naive or random)
    warm = FALSE
    clean.warm.start(warm.start)
    if (!is.null(warm.start)) {
      #must have u,d and v components
      if (!all(match(c("u", "d", "v", "ub", "db", "vb"), names(warm.start), 0) >
               0))
        stop("warm.start is missing some or all of the components.")
      warm = TRUE
      # prepare U, Dsq, V
      D = warm.start$d
      JD = sum(D > min_eigv)
      if (JD >= J) {
        U = warm.start$u[, seq(J), drop = FALSE]
        V = warm.start$v[, seq(J), drop = FALSE]
        Dsq = D[seq(J)]
      } else{
        # upscale
        Ja = J - JD
        Dsq = c(D, rep(D[JD], Ja))
        U = warm.start$u
        Ua = matrix(rnorm(n * Ja), n, Ja)
        Ua = Ua - U %*% (t(U) %*% Ua)
        Ua = tryCatch(
          fast.svd(Ua, trim = FALSE)$u,
          error = function(e)
            svd(Ua)$u
        )
        U = cbind(U, Ua)
        V = cbind(warm.start$v, matrix(0, m, Ja))
      }
      #----------------
      # prepare Q, and R
      Db = warm.start$db
      rD = sum(Db > 0)
      if (rD >= r) {
        Ub = warm.start$ub[, seq(r), drop = FALSE]
        Vb = warm.start$vb[, seq(r), drop = FALSE]
        Db = Db[seq(r)]
      } else{
        ra = r - rD
        Db = c(Db, rep(Db[rD], ra))
        Ub = warm.start$ub
        Uba = matrix(rnorm(k * ra), k, ra)
        Uba = Uba - Ub %*% (t(Ub) %*% Uba)
        Uba = tryCatch(
          fast.svd(Uba, trim = FALSE)$u,
          error = function(e)
            svd(Uba)$u
        )
        Ub = cbind(Ub, Uba)
        Vb = cbind(warm.start$vb, matrix(0, m, ra))
      }
      
      
    } else{
      # initialize. Warm start is not provided
      stopifnot(init %in% c("random", "naive"))
      if (init == "random") {
        stop("Not implemented yet.")
      } else if (init == "naive") {
        Y_naive = as.matrix(y)
        Y_naive = naive_MC(Y_naive)
        svdH = reduced_hat_decomp.H(X)
        Xbeta <-  svdH$u %*% (svdH$v  %*% Y_naive)
        M <- as.matrix(Y_naive - Xbeta)
        M <-  tryCatch(
          propack.svd(M, J),
          error = function(e) {
            message(paste("Naive Init:", e))
            svd_trunc_simple(M, J)
          }
        )
        U = M$u
        V = M$v
        Dsq = pmax(M$d, min_eigv)
        #----------------------
        # initialization for beta = X^-1 Y
        # comment for later: shouldn't be X^-1 H Y??
        beta = as.matrix(ginv(X) %*% Xbeta)
        QRsvd = svd_trunc_simple(beta, r)
        Ub <- as.matrix(QRsvd$u, k, r)
        Db <- sqrt(QRsvd$d)
        Vb <- as.matrix(QRsvd$v, m, r)
        Y_naive <- Xbeta <- M <- NULL
        #---------------------------------------------------------------
      }
    }
    Q = UD(Ub, Db)
    R = UD(Vb, Db)
    XtX = t(X) %*% X
    #----------------------------------------
    yobs <- y@x # y will hold the model residuals
    ratio <- 1
    iter <- 0
    if (!is.null(r) && r == 0) {
      beta <- matrix(0, k, m)
      xbeta.obs <- rep(0, length(y@x))
    }
    #----------------------------------------
    while ((ratio > thresh) & (iter < maxit)) {
      iter <- iter + 1
      U.old = U
      V.old = V
      Dsq.old = Dsq
      Ub.old = Ub
      Vb.old = Vb
      Db.old = Db
      #----------------------------------------------
      # part 0: Update y (training residulas)
      # prereq: U, Dsq, V, Q, R, yobs
      # updates: y; VDsq; XQ
      VDsq = t(Dsq * t(V))
      XQ = as.matrix(X %*% Q)
      M_obs = suvC(U, VDsq, irow, pcol)
      if (iter == 1) {
        xbeta.obs <- suvC(XQ, R, irow, pcol)
        y@x = yobs - M_obs - xbeta.obs
      }
      #----------------------------------------------
      # part 1: Update R
      # prereq: Q, R, XtX, lambda.beta, y, X, r
      # updates: Q, R
      part1 = as.matrix(t(Q) %*% XtX %*% Q + diag(lambda.beta, r, r))
      part2 =  as.matrix(t(XQ) %*% y + (t(Q) %*% (XtX %*% Q)) %*% t(R))
      RD = t(as.matrix(ginv(part1) %*% part2) * Db)
      
      Rsvd = fast.svd(RD, trim = FALSE)
      Ub = Ub %*% Rsvd$v
      Db <- sqrt(Rsvd$d)
      Vb = Rsvd$u
      Q = UD(Ub , Db)
      R = UD(Vb , Db)
      #-----------------------------------------------------------------
      # part 2: update Q
      # prereq: Q, R, X, lambda.beta, XtX, y, Rsvd
      # updates: Q, R
      
      if(Qtype == 0){
        
      part1 <- as.vector(t(X) %*% y %*% R + XtX %*% UD(Q, Db^2))
      part2 <- kronecker(diag(Db^2, r, r), XtX) + diag(lambda.beta, k*r,k*r)
      Q <- matrix(solve(part2) %*% part1, k, r)
      qiter = 0
      }else{
        
      
      # 
      # 
      Qold = Q + 10
      qiter = 0
      Dbsq = Db ^ 2
      # newton!!
      if (Qtype == 1) {
        hessian <- matrix(0, k, r)
        for (ih in 1:k)
          hessian[ih,] <- sum(XtX[ih,]) * Dbsq
        hessian = hessian + diag(lambda.beta, k, r)
        partial_gradient = -t(X) %*% y %*% R - XtX %*% UD(Q, Dbsq)

        while (sqrt(sum((Q - Qold) ^ 2)) > 1e-3 &
               qiter < qiter.max) {
          Qold <- Q
          gradient = partial_gradient +  XtX %*% UD(Q, Dbsq) + lambda.beta * Q
          Q <- Q -  gradient / hessian
          qiter <- qiter + 1
        }
      } else if (Qtype == 2) {
        # proximal
        eta = sum(Dbsq^2) * learning_rate
        #print(eta)
        #print(((1 + lambda.beta * eta)^(-1)))
        partial_grad = t(X) %*% y %*% R + XtX %*% UD(Q, Dbsq)
        while (sqrt(sum((Q - Qold) ^ 2)) > 1e-3 &
               qiter < qiter.max) {
          Qold <- Q
          #print((XtX %*% UD(Q-Qold, Dbsq))[1:3,1:3])
          # print(Q[1:3,1:3])
          Q = ((1 + lambda.beta * eta)^(-1)) *
            (Q + eta * (partial_grad - XtX %*% UD(Q, Dbsq)))
          qiter <- qiter + 1
        }
      }
      #print(qiter)
      }
      
      Qsvd = fast.svd(UD(Q, Db), trim = FALSE)
      Ub = Qsvd$u
      Db <- sqrt(Qsvd$d)
      Vb = Vb %*% Qsvd$v
      Q = UD(Ub,Db)
      R = UD(Vb,Db)
      #-------------------------------------------------------------------
      # part extra: re-update y
      xbeta.obs <-
        suvC(as.matrix(X %*% Q), as.matrix((R)), irow, pcol)
      y@x = yobs - M_obs - xbeta.obs
      #--------------------------------------------------------------------
      # print("hi2")
      ##--------------------------------------------
      # part 3: Update B
      # prereq: U, VDsq, y, Dsq, lambda.M, L.b
      # updates U, Dsq, V
      B = t(U) %*% y + t(VDsq)
      if (laplace.b)
        B = B - t(V)  %*% L.b
      B = as.matrix(t((B) * (Dsq / (Dsq + lambda.M))))
      Bsvd = tryCatch({
        fast.svd(B, trim = FALSE)
      }, error = function(e) {
        message(paste("Loop/B:", e))
        svd(B)
      })
      V = Bsvd$u
      Dsq = pmax(Bsvd$d, min_eigv)
      U = U %*% (Bsvd$v)
      #-------------------------------------------------------------
      # part 4: Update A
      # prereq: U, D, VDsq, y, lambda.M, L.a
      # updates U, Dsq, V
      A = (y %*% V) + t(Dsq * t(U))
      if (laplace.a)
        A = A - L.a %*% U
      A = as.matrix(t(t(A) * (Dsq / (Dsq + lambda.M))))
      Asvd = tryCatch({
        #svd(A)
        fast.svd(A, trim = FALSE)
      }, error = function(e) {
        message(paste("Loop/A:", e))
        svd(A)
      })
      U = Asvd$u
      Dsq = pmax(Asvd$d, min_eigv)
      V = V %*% (Asvd$v)
      #------------------------------------------------------------------------------
      ratio =  Frob(U.old, Dsq.old, V.old, U, Dsq, V)
      ratio = ratio + Frob(Ub.old, Db.old, Vb.old, Ub, Db, Vb)
      ratio = ratio / 2
      #------------------------------------------------------------------------------
      if (trace.it) {
        obj = (.5 * sum(y@x ^ 2) +
                 lambda.M * sum(Dsq) +
                 lambda.beta * sum(Db ^ 2)) / nz
        cat(iter,
            ":",
            "obj",
            format(round(obj, 5)),
            "ratio",
            ratio,
            "qiter",
            qiter,
            "\n")
      }
      #-----------------------------------------------------------------------------------
      
    }
    if (iter == maxit &
        trace.it)
      warning(
        paste(
          "Convergence not achieved by",
          maxit,
          "iterations. Consider increasing the number of iterations."
        )
      )
    # one final fit for one of the parameters (A) has proved to improve the performance significantly.
    if (final.svd) {
      #---- update A
      A = (y %*% V) + t(Dsq * t(U))
      if (laplace.a)
        A = A - L.a %*% U
      A = as.matrix(t(t(A) * (Dsq / (Dsq + lambda.M))))
      Asvd = tryCatch({
        fast.svd(A, trim = FALSE)
      }, error = function(e) {
        message(paste("Final/A:", e))
        svd(A)
      })
      U = Asvd$u
      V = V %*% Asvd$v
      # is this good?
      Dsq = pmax(Asvd$d - lambda.M, min_eigv)
      #---------------------------------------------
      if (trace.it) {
        M_obs = suvC(t(Dsq * t(U)), V, irow, pcol)
        y@x = yobs - M_obs - xbeta.obs
        obj = (.5 * sum(y@x ^ 2) + lambda.M * sum(Dsq)) / nz
        cat("final SVD:", "obj", format(round(obj, 5)), "\n")
        cat("Number of Iterations for covergence is ", iter, "\n")
      }
    }
    
    #-------------------------------------------------------
    # trim in case we reduce the rank of M to be smaller than J.
    J = min(sum(Dsq > min_eigv) + 1, J)
    J = min(J, length(Dsq))
    r = min(sum(Db > 0) + 1, r)
    r = min(r, length(Db))
    
    out = list(
      u = U[, seq(J), drop = FALSE],
      d = Dsq[seq(J)],
      v = V[, seq(J), drop = FALSE],
      ub = Ub[, seq(r), drop = FALSE],
      db = Db[seq(r)],
      vb = Vb[, seq(r), drop = FALSE],
      #---------------------
      # hyperparameters
      lambda.M = lambda.M,
      lambda.beta = lambda.beta,
      J = J,
      r = r,
      lambda.a = lambda.a,
      lambda.b = lambda.b,
      #--------------------------
      # convergence
      n_iter = iter
    )
    out
  }

