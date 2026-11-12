# Independent base-R dense REML and direct GLS score calculations.
# Run from repository root. Does not call JLinAlg or GPU_eQTL.
options(digits=17)
out <- "test/resources/related-set-reference"
dir.create(out, recursive=TRUE, showWarnings=FALSE)
d <- read.csv("test/resources/mixed-reference/data.csv")
n <- nrow(d); X <- cbind(1,d$x); A <- diag(n)
for(f in seq(1,n,4)) for(child in c(f+2,f+3))
  for(j in seq_len(child-1)) A[child,j] <- A[j,child] <- (A[f,j]+A[f+1,j])/2
y <- (d$y-min(d$y))/(max(d$y)-min(d$y)) # observed endpoints 0 and 1
G <- cbind((seq_len(n)%%7==0)*1,(seq_len(n)%%5==0)*1,(seq_len(n)%%4==0)*1)
G[c(6,14),1] <- NA; G[24,2] <- NA
write.csv(data.frame(ID=d$id,x=d$x,y=y,g1=G[,1],g2=G[,2],g3=G[,3]),file.path(out,"data.csv"),row.names=FALSE,quote=FALSE,na="NA")
write.table(data.frame(ID=d$id,A,check.names=FALSE),file.path(out,"grm.tsv"),sep="\t",row.names=FALSE,quote=FALSE,col.names=c("ID",d$id))
fit <- function(y,X,K) {
  bases <- list(K,diag(length(y))); n <- length(y)
  evaluate <- function(logvar,details=FALSE) {
    v <- exp(logvar); V <- v[1]*K+v[2]*diag(n); L <- chol(V); vi <- chol2inv(L)
    I <- crossprod(X,vi%*%X); P <- vi-vi%*%X%*%solve(I)%*%t(X)%*%vi
    py <- drop(P%*%y)
    loss <- .5*((n-ncol(X))*log(2*pi)+2*sum(log(diag(L)))+as.numeric(determinant(I,logarithm=TRUE)$modulus)+sum(y*py))
    gradient <- -.5*v*vapply(bases,function(k) sum(py*(k%*%py))-sum(diag(P%*%k)),0.)
    if(details) list(P=P,var=v,loglik=-loss,gradient=gradient) else loss
  }
  opt <- optim(rep(log(var(y)/2),2),evaluate,gr=function(v)evaluate(v,TRUE)$gradient,
    method="L-BFGS-B",lower=log(var(y)*1e-10),upper=log(var(y)*1e10),
    control=list(factr=1,pgtol=1e-10,maxit=3000))
  stopifnot(opt$convergence==0); evaluate(opt$par,TRUE)
}
expected <- data.frame()
for(case in c("complete","pattern","same-mask-new-trait")) {
  mask <- if(case=="complete") seq_len(n) else setdiff(seq_len(n),c(3,7,11,20,29))
  yy <- if(case=="same-mask-new-trait") y[mask]^2 else y[mask]
  gg <- G[mask,,drop=FALSE]
  for(j in 1:3) gg[is.na(gg[,j]),j] <- mean(gg[,j],na.rm=TRUE)
  gg <- sweep(gg,2,c(1,.7,1.3),`*`)
  for(related in c(FALSE,TRUE)) {
    if(related) {
      null <- fit(yy,X[mask,],A[mask,mask]); P <- null$P
      U <- drop(crossprod(gg,P%*%yy)); C <- crossprod(gg,P%*%gg)
      variances <- null$var; ll <- null$loglik
    } else {
      P <- diag(length(mask))-X[mask,]%*%solve(crossprod(X[mask,]))%*%t(X[mask,])
      sigma <- sqrt(sum(yy*(P%*%yy))/(length(mask)-ncol(X)))
      U <- drop(crossprod(gg,P%*%yy))/sigma; C <- crossprod(gg,P%*%gg)
      variances <- c(NA,NA); ll <- NA
    }
    expected <- rbind(expected,data.frame(case=case,related=related,N=length(mask),DF=length(mask)-ncol(X),
      u1=U[1],u2=U[2],u3=U[3],c11=C[1,1],c12=C[1,2],c13=C[1,3],c22=C[2,2],c23=C[2,3],c33=C[3,3],
      burden=sum(U)^2/sum(C),skat=sum(U^2),genetic=variances[1],residual=variances[2],loglik=ll))
  }
}
write.csv(expected,file.path(out,"scores.csv"),row.names=FALSE,quote=FALSE)
writeLines(trimws(c("Base R optim/chol/solve oracle; exact mask subset before REML; mean fill after subset.",
  capture.output(sessionInfo())),which="right"),file.path(out,"session-info.txt"))
