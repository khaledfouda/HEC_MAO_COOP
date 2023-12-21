SVT_Approx <- function(A, R, lambda, maxIter=100, tol=1e-6) {
   
   
   # the following part computes Algorithm 1: Compute Q
   Q <- qr(A %*% R, LAPACK = TRUE)$Q
   
   err <- numeric(maxIter)
   for (i in seq_len(maxIter)) {
      iQ <- qr(A %*% t(A %*% Q), LAPACK = TRUE)$Q
      
      err[i] <- norm(iQ[,1] - Q[,1], type = "2")
      Q <- iQ
      
      if(err[i] < tol) {
         break
      }
   }
   #------------
   # the next part computes Algorithm 2:  Stage 2
   svdQA = svd(t(Q)%*%A)
   j = sum( (svdQA$d-lambda) > 0 )
   return(list(u = (Q%*%svdQA$u)[,1:j], d=svdQA$d[1:j], v=svdQA$v[,1:j]))
}
