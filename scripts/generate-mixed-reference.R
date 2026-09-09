# Independent Gaussian REML reference using base R Cholesky and optim, plus nlme.
# No JLinAlg or GPU eQTL code is used to generate the expected values.
options(digits=17)
out <- "test/resources/mixed-reference"
dir.create(out, recursive=TRUE, showWarnings=FALSE)
n <- 48L
id <- paste0("s",1:n)
sire <- dam <- rep("0",n)
A <- diag(n)
for(family in 0:11) {
  f <- 4*family+1; m <- f+1
  for(child in c(f+2,f+3)) {
    sire[child] <- id[f]; dam[child] <- id[m]
    for(j in 1:(child-1)) A[child,j] <- A[j,child] <- .5*(A[f,j]+A[m,j])
  }
}
group <- rep(1:8,each=6)
K <- outer(group,group,"==")*1
x <- sin((1:n)*.73)
g <- rep(c(0,1,2,1,0,2),8)+.1*cos(1:n)
X <- cbind(1,x)
V <- 2*A + 3*K + .4*diag(n)
set.seed(73091)
y <- 1.3 + .4*x + .35*g + drop(t(chol(V))%*%rnorm(n))
write.csv(data.frame(id,group,x,g,y),file.path(out,"data.csv"),row.names=FALSE,quote=FALSE)
write.table(data.frame(ID=id,SIRE=sire,DAM=dam),file.path(out,"pedigree.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
reml <- function(X, bases) {
  evaluate <- function(logvar, details=FALSE) {
    variance <- exp(logvar)
    V <- Reduce(`+`,Map(`*`,bases,variance))
    cholV <- chol(V); vi <- chol2inv(cholV)
    information <- crossprod(X,vi%*%X)
    bc <- solve(information)
    beta <- drop(bc%*%crossprod(X,vi%*%y))
    resid <- y-drop(X%*%beta)
    ll <- -.5*((n-ncol(X))*log(2*pi)+2*sum(log(diag(cholV)))+
       as.numeric(determinant(information,logarithm=TRUE)$modulus)+sum(resid*(vi%*%resid)))
    if(!details) return(-ll)
    P <- vi-vi%*%X%*%bc%*%t(X)%*%vi
    py <- drop(P%*%y)
    gradient <- -.5*variance*vapply(bases,function(K) sum(py*(K%*%py))-sum(diag(P%*%K)),0.)
    list(var=variance,beta=beta,bc=bc,loglik=ll,P=P,gradient=gradient)
  }
  opt <- optim(rep(log(var(y)/length(bases)),length(bases)),evaluate,
    gr=function(v)evaluate(v,TRUE)$gradient,method="L-BFGS-B",lower=log(1e-10),upper=log(1e10),
    control=list(factr=1e3,pgtol=1e-8,maxit=3000))
  stopifnot(opt$convergence==0)
  evaluate(opt$par,TRUE)
}
result <- data.frame()
for (name in c("group","pedigree","combined")) {
  bases <- switch(name,group=list(K,diag(n)),pedigree=list(A,diag(n)),combined=list(A,K,diag(n)))
  null <- reml(X,bases); full <- reml(cbind(X,g),bases)
  U <- sum(g*(null$P%*%y)); information <- sum(g*(null$P%*%g))
  score <- U^2/information; beta <- U/information; se <- 1/sqrt(information)
  result <- rbind(result,data.frame(model=name,mode="score",beta,se,statistic=score,df=1,
    log10p=pchisq(score,1,lower.tail=FALSE,log.p=TRUE)/log(10),loglik=null$loglik,
    variances=paste(null$var,collapse=";")))
  beta <- tail(full$beta,1); se <- sqrt(tail(diag(full$bc),1)); stat <- beta/se; df <- n-ncol(X)-1
  result <- rbind(result,data.frame(model=name,mode="reml",beta,se,statistic=stat,df,
    log10p=(log(2)+pt(abs(stat),df,lower.tail=FALSE,log.p=TRUE))/log(10),loglik=full$loglik,
    variances=paste(full$var,collapse=";")))
  if(name=="group") {
    oracle <- nlme::lme(y~x+g,random=~1|group,method="REML",control=nlme::lmeControl(msTol=1e-10))
    stopifnot(max(abs(nlme::fixef(oracle)-full$beta))<1e-6, abs(as.numeric(logLik(oracle))-full$loglik)<1e-7)
  }
}
write.csv(result,file.path(out,"association.csv"),row.names=FALSE,quote=FALSE)
writeLines(sub("[[:blank:]]+$", "", capture.output(sessionInfo())),file.path(out,"session-info.txt"))
