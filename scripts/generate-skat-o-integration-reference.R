# Independent reference: R stats::integrate/uniroot with direct chi-square tails.
options(digits=17)
.libPaths(c("target/rare-r-lib",.libPaths()))
out <- "test/resources/related-set-reference"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
mixture2 <- function(q,lambda) {
  if(q==0) return(1)
  peak <- exp(-q/(2*max(lambda)))
  peak*integrate(function(t) exp(-q/2*(1/(lambda[1]*cos(t)^2+lambda[2]*sin(t)^2)-1/max(lambda))),
    0,pi/2,rel.tol=1e-11,abs.tol=1e-13,subdivisions=4000)$value*2/pi
}
critical <- function(p,lambda) {
  if(length(lambda)==1) return(lambda*qchisq(p,1,lower.tail=FALSE))
  uniroot(function(q) log(mixture2(q,lambda))-log(p),c(0,max(lambda)*qchisq(p,2,lower.tail=FALSE)),tol=1e-11)$root
}
records <- data.frame()
for(cross in c(.25,-.4)) for(a in c(3,6,9,12)) {
  C <- matrix(c(1,cross,cross,.7),2); U <- a*c(.8,-.4); L <- t(chol(C))
  rho <- c(0,.25,.5,.75,1)
  matrices <- lapply(rho,function(r) t(L)%*%((1-r)*diag(2)+r*matrix(1,2,2))%*%L)
  spectra <- lapply(matrices,function(M) { x<-eigen(M,symmetric=TRUE,only.values=TRUE)$values; x[x>1e-12] })
  q <- (1-rho)*sum(U^2)+rho*sum(U)^2
  p <- mapply(function(q,l) if(length(l)==1) pchisq(q/l,1,lower.tail=FALSE) else mixture2(q,l),q,spectra)
  minimum <- min(p); crit <- vapply(spectra,function(l)critical(minimum,l),0.)
  # Condition on the first standard normal. Intersection of non-rejection
  # intervals in the second normal is [max(lower),min(upper)]. Integrate its
  # complement using direct normal tails, never 1-nearly-one.
  conditional <- function(z) vapply(z,function(z1) {
    low <- -Inf; high <- Inf
    for(i in seq_along(rho)) {
      M <- matrices[[i]]; aa <- M[2,2]; bb <- 2*M[1,2]*z1; cc <- M[1,1]*z1^2-crit[i]
      discriminant <- bb^2-4*aa*cc
      if(discriminant<=0) return(dnorm(z1)/minimum)
      roots <- c((-bb-sqrt(discriminant))/(2*aa),(-bb+sqrt(discriminant))/(2*aa))
      low <- max(low,roots[1]); high <- min(high,roots[2])
    }
    if(low>=high) return(dnorm(z1)/minimum)
    (pnorm(low)+pnorm(high,lower.tail=FALSE))*dnorm(z1)/minimum
  },0.)
  # Partition to resolve all changes of the limiting component.
  value <- sum(vapply(seq(-14,13.5,.5),function(z)integrate(conditional,z,z+.5,
    rel.tol=1e-9,abs.tol=1e-11,subdivisions=4000)$value,0.))*minimum
  records <- rbind(records,data.frame(cross=cross,amplitude=a,minimum=minimum,adjusted=value))
}
write.csv(records,file.path(out,"skat-o-rank-two.csv"),row.names=FALSE,quote=FALSE)
# Orthogonal equal-variance scores, rho={0,1}. Exact joint event is
# {chi-square(m)>c0 OR Z_burden^2>c1}; use independent residual chi-square(m-1).
records <- data.frame()
for(m in c(3,8,20)) for(p in c(1e-4,1e-6,1e-8,1e-10,1e-12)) {
  c0 <- qchisq(p,m,lower.tail=FALSE); c1 <- qchisq(p,1,lower.tail=FALSE)
  f <- function(z) 2*dnorm(z)*pchisq(pmax(0,c0-z*z),m-1,lower.tail=FALSE)/p
  adjusted <- p+integrate(f,0,sqrt(c1),rel.tol=1e-10,abs.tol=1e-12,subdivisions=4000)$value*p
  # Scores perpendicular to burden, so only SKAT determines observed minimum p.
  u <- c(sqrt(c0/2),-sqrt(c0/2),rep(0,m-2))
  records <- rbind(records,data.frame(rank=m,minimum=p,score=sqrt(c0/2),adjusted=adjusted))
}
write.csv(records,file.path(out,"skat-o-isotropic.csv"),row.names=FALSE,quote=FALSE)
writeLines(trimws(c("stats::integrate (QUADPACK) direct conditional normal / chi-square tails",
  capture.output(sessionInfo())),which="right"),file.path(out,"skat-o-session-info.txt"))
