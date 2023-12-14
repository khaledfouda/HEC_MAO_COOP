Frob=function(Uold,Dsqold,Vold,U,Dsq,V){
   denom=sum(Dsqold^2)
   utu=Dsq* (t(U)%*%Uold)
   vtv=Dsqold* (t(Vold)%*%V)
   uvprod= sum(diag(utu%*%vtv))
   num=denom+sum(Dsq^2) -2*uvprod
   num/max(denom,1e-9)
}

clean.warm.start=function(a){
   if(is.null(a))return(NULL)
   d=a$d
   if(is.null(d))return(NULL)
   if(any(d>0)){
      if(length(d)==1){
         a$u=matrix(a$u,ncol=1)
         a$v=matrix(a$v,ncol=1)
      }
      a
   }
   else NULL
}

simpute.als.cov <-
   function (Y, X, beta_partial, J = 2, thresh = 1e-05,lambda=0,maxit=100,trace.it=FALSE,warm.start=NULL){
   
   # are you scaling???
   
   n1 <- dim(Y)[1]
   n2 <- dim(Y)[2]
   m1 <- dim(X)[2]
   ynas <- Y == 0
   #ynas <- is.na(Y)
   nz=n1*n2-sum(ynas)
   
   # The following two lines are as shown in (c) and (d)
   yfill <- Y #as(Y, "sparseMatrix")
   #yfill[ynas] <- 0
   
   
   if(!is.null(warm.start)){
      ###must have u,d and v components
      warm.start = clean.warm.start(warm.start)
      #-------------------------------------------------------------------
      if(!all(match(c("u","d","v", "beta.estim"),names(warm.start),0)>0))stop("warm.start does not have components u, d and v")
      warm=TRUE
      D=warm.start$d
      JD=sum(D>0)
      Xbeta <- X %*% warm.start$beta.estim
      if(JD >= J){
         U=warm.start$u[,seq(J),drop=FALSE]
         V=warm.start$v[,seq(J),drop=FALSE]
         Dsq=D[seq(J)]
      }
      else{
         Dsq=c(D,rep(D[JD],J-JD))
         Ja=J-JD
         U=warm.start$u
         Ua=matrix(rnorm(n1*Ja),n1,Ja)
         Ua=Ua-U%*% (t(U)%*%Ua)
         Ua=svd(Ua)$u
         U=cbind(U,Ua)
         V=cbind(warm.start$v,matrix(0,n2,Ja))
      }
      yfill[ynas]=(U%*%(Dsq*t(V)))[ynas] + Xbeta[ynas]  
   }
   else
   {
      V=matrix(0,n2,J)
      U=matrix(rnorm(n1*J),n1,J)
      U=svd(U)$u
      Dsq=rep(1,J)# we call it Dsq because A=UD and B=VD and AB'=U Dsq V^T
      yfill[ynas]=0
   }
      
   beta.estim <- beta_partial %*% yfill
   Xbeta <- X %*% beta.estim
   yplus <- yfill - Xbeta 
   ratio <- 1
   iter <- 0
   while ((ratio > thresh)&(iter<maxit)) {
      
      iter <- iter + 1
      U.old=U
      V.old=V
      Dsq.old=Dsq
      ## U step
      B=t(U)%*%yplus
      if(lambda>0)B=B*(Dsq/(Dsq+lambda))
      Bsvd=svd(t(B))
      V=Bsvd$u
      Dsq=(Bsvd$d)
      U=U%*%Bsvd$v
      yhat=U %*%(Dsq*t(V)) # yhat = AB - Xbeta
      yfill[ynas]=yhat[ynas] + Xbeta[ynas]
      # new beta estimates
      beta.estim <- beta_partial %*% yfill
      Xbeta <- X %*% beta.estim
      yplus <- yfill - Xbeta # updated_Y - Xbeta
      ###The next line we could have done later; this is to match with sparse version
      if(trace.it) obj=(.5*sum( (yplus-yhat)[!ynas]^2)+lambda*sum(Dsq))/nz
      ## V step
      A=t(yplus%*%V)
      if(lambda>0)A=A*(Dsq/(Dsq+lambda))
      Asvd=svd(t(A))
      U=Asvd$u
      Dsq=Asvd$d
      V=V %*% Asvd$v
      yhat=U %*%(Dsq*t(V)) # yhat = AB - Xbeta
      yfill[ynas]=yhat[ynas] + Xbeta[ynas]
      # new beta estimates
      beta.estim <- beta_partial %*% yfill
      Xbeta <- X %*% beta.estim
      yplus <- yfill - Xbeta # yplus = updated_Y (yfill) - Xbeta
      ratio=Frob(U.old,Dsq.old,V.old,U,Dsq,V)
      if(trace.it) cat(iter, ":", "obj",format(round(obj,5)),"ratio", ratio, "\n")
      ##########################################################
   }
   if(iter==maxit) warning(paste("Convergence not achieved by",maxit,"iterations"))
   
   U=yplus%*%V
   sU=svd(U)
   U=sU$u
   Dsq=sU$d
   V=V%*%sU$v
   Dsq=pmax(Dsq-lambda,0)
   if(trace.it){
      yhat=U %*%(Dsq*t(V))
      obj=(.5*sum( (yplus-yhat)[!ynas]^2)+lambda*sum(Dsq))/nz
      cat("final SVD:", "obj",format(round(obj,5)),"\n")
   }
   J=min(sum(Dsq>0)+1,J)
   out=list(u=U[,seq(J)],d=Dsq[seq(J)],v=V[,seq(J)],J=J, lambda=lambda, beta.estim=beta.estim)
   out
}
######################



# 
# set.seed(101)
# n=200
# p=100
# J=50
# np=n*p
# missfrac=0.3
# x=matrix(rnorm(n*J),n,J)%*%matrix(rnorm(J*p),J,p)+matrix(rnorm(np),n,p)/5
# ix=seq(np)
# imiss=sample(ix,np*missfrac,replace=FALSE)
# xna=x
# xna[imiss]=NA
# xnaC=as(xna,"Incomplete")
# ### here we do it a different way to demonstrate Incomplete
# ### In practise the observed values are stored in this market-matrix format.
# i = row(xna)[-imiss]
# j = col(xna)[-imiss]
# xnaC=Incomplete(i,j,x=x[-imiss])