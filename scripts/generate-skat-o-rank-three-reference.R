# Established CompQuadForm::farebrother mixture tails + cubature::hcubature
# integration of the joint Gaussian radial event; no GPU_eQTL/JLinAlg code.
.libPaths(c("target/rare-r-lib",.libPaths())); options(digits=17)
lambda <- c(.8,1,1.2); rho <- c(0,.25,.5,.75,1); b <- sqrt(lambda)
matrices <- lapply(rho,function(r)(1-r)*diag(lambda)+r*tcrossprod(b))
spectra <- lapply(matrices,function(M){v<-eigen(M,symmetric=TRUE)$values; v[v>1e-12]})
tail <- function(q,l) {
  if(q==0) return(1)
  if(length(l)==1) return(pchisq(q/l,1,lower.tail=FALSE))
  ans <- CompQuadForm::farebrother(q,l,eps=1e-15,maxit=100000)
  stopifnot(ans$ifault==0); ans$Qq
}
quantile <- function(p,l) {
  if(length(l)==1) return(l*qchisq(p,1,lower.tail=FALSE))
  uniroot(function(q)tail(q,l)-p,c(0,1.01*max(l)*qchisq(p,length(l),lower.tail=FALSE)),tol=1e-10)$root
}
result <- data.frame()
for(p in c(1e-4,1e-6,1e-8,1e-10)) {
  crit <- vapply(spectra,function(l)quantile(p,l),0.)
  fn <- function(x) {
    z <- x[1]; phi <- x[2]; d <- c(sqrt(1-z*z)*cos(phi),sqrt(1-z*z)*sin(phi),z)
    qs <- sum(lambda*d*d); qb <- sum(b*d)^2
    boundary <- min(crit/((1-rho)*qs+rho*qb))
    pchisq(boundary,3,lower.tail=FALSE)/(4*pi*p)
  }
  integrated <- cubature::hcubature(fn,c(-1,0),c(1,2*pi),tol=1e-6,maxEval=4000000)
  stopifnot(integrated$returnCode==0,integrated$error/integrated$integral < 2e-6)
  result <- rbind(result,data.frame(minimum=p,score=sqrt(crit[1]/2),adjusted=p*integrated$integral,
    estimated_error=p*integrated$error,evaluations=integrated$functionEvaluations))
}
write.csv(result,"test/resources/related-set-reference/skat-o-rank-three.csv",row.names=FALSE,quote=FALSE)
writeLines(trimws(c(paste("CompQuadForm",packageVersion("CompQuadForm")),paste("cubature",packageVersion("cubature")),
  capture.output(sessionInfo())),which="right"),"test/resources/related-set-reference/rank-three-session-info.txt")
