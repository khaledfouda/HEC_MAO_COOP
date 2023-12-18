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

simpute.als.cov.v2 <-
   function (Y, X, beta_partial, J = 2, thresh = 1e-05,lambda=0,maxit=100,trace.it=FALSE,warm.start=NULL){
   
   # are you scaling???
   
   n1 <- dim(Y)[1]
   n2 <- dim(Y)[2]
   m1 <- dim(X)[2]
   ymiss <- Y == 0
   yobs <- ! ymiss
   
   nz=n1*n2-sum(ymiss)
   
   # The following two lines are as shown in (c) and (d)
   #yfill <- Y #as(Y, "sparseMatrix")
   #yfill[ymiss] <- 0
   
   Sp <- Y
   
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
      yfill[ymiss]=(U%*%(Dsq*t(V)))[ymiss] #+ Xbeta[ymiss]  
   }
   else
   {
      V=matrix(0,n2,J)
      U=matrix(rnorm(n1*J),n1,J)
      U=svd(U)$u
      Dsq=rep(1,J)# we call it Dsq because A=UD and B=VD and AB'=U Dsq V^T
      #yfill[ymiss]=0
   }
   BD = V * Dsq
   AD <- U * Dsq
   
   AB = (U * Dsq) %*% t(V)
   
   yfill <- Y
   beta.estim <- beta_partial %*% yfill
   Xbeta <- X %*% beta.estim
   Sp[yobs] = Y[yobs] - AB[yobs] - Xbeta[yobs]
   #beta.estim <- beta_partial %*% yfill
   #Xbeta <- X %*% beta.estim
   #yplus <- yfill #- Xbeta 
   #Sp[yobs] = Y[yobs] - AB[yobs] #- Xbeta[yobs]
    #  X.star = Sp + AB
   ratio <- 1
   iter <- 0
   while ((ratio > thresh)&(iter<maxit)) {
      
      iter <- iter + 1
      U.old=U
      V.old=V
      Dsq.old=Dsq
      #-------------------------------------------------------
      # Need:  Sp, Dsq, U, BD
      # Estimate new BD
      # this is the part P = solve( D%*%D + lambda I) %*% D %*% D
      # but since D will always be diagonal, the following is faster
      # We also keeep Dsq as a vector for optimization
      # We will actually estimate BD instead of B because 
      #   1. We will always need B^2 not B and 2. we need BD not B
      #BD = t(P %*% t(U) %*% Sp + P %*% t(BD))
      #this one matches theirs for some reason
      #BD = t((t(U)%*%X.star)*P)
      #BD = t((t(U)%*%Sp)*P + t(U) %*% U %*% t(BD)*P)
      #BD = t(Sp) %*% U %*% P + BD %*% P
      P = ( Dsq / (Dsq + lambda))
      BD = t((t(U)%*%Sp)*P + (  t(V) * (Dsq *P)) + (t(U)%*% Xbeta)*P )  
      # Now we do SVD decomposition to compute new V, and D values
      svd.BD = svd(BD)
      V = svd.BD$u
      Dsq <- svd.BD$d
      BD = V * Dsq 
      U = U %*% svd.BD$v
      #------------------------------------------------------------
      # Update Sp 
      AB = (U) %*% (Dsq * t(V)) 
      #-----------------------------------------------
      # Estimate Beta and update Sp
      yfill[ymiss] <- AB[ymiss] + Xbeta[ymiss]
      beta.estim <- beta_partial %*% yfill
      Xbeta <- X %*% beta.estim
      Sp[yobs] = Y[yobs] - AB[yobs] - Xbeta[yobs]
      #X.star = Sp + AB
      #-------------------------------------------------------------
      # Estimate new AD
      # need: Sp, V, Dsq, AD
      P = ( Dsq / (Dsq + lambda)) # note P is a diagonal
      #AD = t(P %*% t(V) %*% t(Sp) + P %*% t(AD))
      # again same as what the did
      #AD = t((t(V) %*% t(X.star)) * P  )
      AD = t( (t(V) %*% t(Sp)) * P + (t(U) * (Dsq*P))  + (t(V)%*% t(Xbeta))*P ) 
      #AD = Sp %*% V %*% P + AD %*% P
      svd.AD = svd(AD)
      U = svd.AD$u
      Dsq <- svd.AD$d
      AD <- U * Dsq
      V <- V %*% svd.AD$v
      #-----------------------------------------------------
      # next we estimate  Sp =  Y - AB for non-missing only!!! like initialize at Y and keep updating.
      # these are the estimates A %*% t(B)
      AB = (U) %*% (Dsq * t(V))
      # Estimate Beta and update Sp
      yfill[ymiss] <- AB[ymiss] + Xbeta[ymiss]
      beta.estim <- beta_partial %*% yfill
      Xbeta <- X %*% beta.estim
      Sp[yobs] = Y[yobs] - AB[yobs] - Xbeta[yobs]
      
      #X.star = Sp + AB
      #------------------------------------------
      # Estimating Covariate effects
      #yfill[ymiss] = AB[ymiss] #+ Xbeta[ymiss]
      #beta.estim <- beta_partial %*% yfill
      #Xbeta <- X %*% beta.estim
      
      #-------------------------------------------------------
      ###The next line we could have done later; this is to match with sparse version
      if(trace.it) obj=(.5*sum( Sp[yobs]^2)+lambda*sum(Dsq))/nz
      #---------------------------------------------------------
      ratio=Frob(U.old,Dsq.old,V.old,U,Dsq,V)
      if(trace.it) cat(iter, ":", "obj",format(round(obj,5)),"ratio", ratio, "\n")
      #---------------------------------------------------------
      ## U step
      # B=t(U)%*%yplus
      # if(lambda>0)B=B*(Dsq/(Dsq+lambda))
      # Bsvd=svd(t(B))
      # V=Bsvd$u
      # Dsq=(Bsvd$d)
      # U=U%*%Bsvd$v
      # yhat=U %*%(Dsq*t(V)) # yhat = AB - Xbeta
      # yfill[ymiss]=yhat[ymiss] + Xbeta[ymiss]
      # # new beta estimates
      # beta.estim <- beta_partial %*% yfill
      # Xbeta <- X %*% beta.estim
      # yplus <- yfill - Xbeta # updated_Y - Xbeta
      # ###The next line we could have done later; this is to match with sparse version
      # if(trace.it) obj=(.5*sum( (yplus-yhat)[!ymiss]^2)+lambda*sum(Dsq))/nz
      # ## V step
      # A=t(yplus%*%V)
      # if(lambda>0)A=A*(Dsq/(Dsq+lambda))
      # Asvd=svd(t(A))
      # U=Asvd$u
      # Dsq=Asvd$d
      # V=V %*% Asvd$v
      # yhat=U %*%(Dsq*t(V)) # yhat = AB - Xbeta
      # yfill[ymiss]=yhat[ymiss] + Xbeta[ymiss]
      # # new beta estimates
      # beta.estim <- beta_partial %*% yfill
      # Xbeta <- X %*% beta.estim
      # yplus <- yfill - Xbeta # yplus = updated_Y (yfill) - Xbeta
      # ratio=Frob(U.old,Dsq.old,V.old,U,Dsq,V)
      # if(trace.it) cat(iter, ":", "obj",format(round(obj,5)),"ratio", ratio, "\n")
      ##########################################################
   }
   if(iter==maxit) warning(paste("Convergence not achieved by",maxit,"iterations"))
   
   #U=yplus%*%V
   sU=svd(U)
   U=sU$u
   Dsq=sU$d
   V=V%*%sU$v
   Dsq=pmax(Dsq-lambda,0)
   if(trace.it){
      AB=U %*%(Dsq*t(V))
      #yfill[ymiss] = AB[ymiss] + Xbeta[ymiss]
      beta.estim <- 0#beta_partial %*% yfill
      #Xbeta <- X %*% beta.estim
      Sp = Y - AB #- Xbeta
      obj=(.5*sum( Sp[yobs]^2)+lambda*sum(Dsq))/nz
      cat("final SVD:", "obj",format(round(obj,5)),"\n")
   }
   J=min(sum(Dsq>0)+1,J)
   out=list(u=U[,seq(J)],d=Dsq[seq(J)],v=V[,seq(J)],J=J, lambda=lambda, beta.estim=beta.estim)
   out
}
######################
system.time({sout2 <- simpute.als.cov.v2(Y_train, gen.dat$X, beta_partial, 15, 1e-3, 30,trace.it = TRUE)})
system.time({sout1 <- simpute.als.cov(Y_train, gen.dat$X, beta_partial, 15, 1e-3, 30,trace.it = TRUE)})

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