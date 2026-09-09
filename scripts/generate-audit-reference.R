# Independent deterministic reference generation. R 4.6.1, CompQuadForm 1.4.4.
dir.create("test/resources/audit-reference", recursive=TRUE, showWarnings=FALSE)
out <- "test/resources/audit-reference"
options(digits=17)
n <- 12
age <- c(1,3,2,6,4,5,9,7,8,12,10,11)
g <- c(0,1,2,1,0,2,0,1,2,0,2,1)
y <- c(2,1,5,3,2,6,3,4,8,4,8,6)
X <- cbind(1,age)
subject <- c(1,2,1,3,2,1,4,3,4,5,5,5)
# Pedigree relationship matrix from two unrelated parents and ten offspring.
A <- diag(n)
for (i in 3:n) {
  for (j in 1:(i-1)) A[i,j] <- A[j,i] <- (A[1,j]+A[2,j])/2
  A[i,i] <- 1 + A[1,2]/2
}
V <- list(fixed=diag(n), repeated=0.35*outer(subject,subject,"==")+0.65*diag(n), pedigree=0.7*A+0.3*diag(n))
write.csv(data.frame(id=paste0("s",1:n),age,g,y,subject),file.path(out,"data.csv"),row.names=FALSE,quote=FALSE)
refs <- data.frame()
for (name in names(V)) {
  v <- V[[name]]
  vi <- solve(v)
  P <- vi - vi %*% X %*% solve(t(X)%*%vi%*%X) %*% t(X)%*%vi
  xx <- drop(t(g)%*%P%*%g); xy <- drop(t(g)%*%P%*%y); yy <- drop(t(y)%*%P%*%y)
  beta <- xy/xx; r2 <- xy*xy/(xx*yy); df <- n-ncol(X)-1
  tt <- beta/sqrt((yy-xy*xy/xx)/df/xx)
  refs <- rbind(refs,data.frame(model=name,beta,r2,t=tt,df,log10p=(log(2)+pt(abs(tt),df,lower.tail=FALSE,log.p=TRUE))/log(10)))
  write.table(data.frame(ID=paste0("s",1:n),v,check.names=FALSE),file.path(out,paste0(name,".tsv")),sep="\t",quote=FALSE,row.names=FALSE,col.names=c("ID",paste0("s",1:n)))
}
write.csv(refs,file.path(out,"gls.csv"),row.names=FALSE,quote=FALSE)
fixedP <- diag(n)-X%*%solve(t(X)%*%X)%*%t(X)
score <- drop(t(g)%*%fixedP%*%y)
information <- drop(t(g)%*%fixedP%*%g)
sigma2 <- drop(t(y)%*%fixedP%*%y)/(n-ncol(X))
write.csv(data.frame(q=score^2/sigma2,lambda=information,p=pchisq(score^2/sigma2/information,1,lower.tail=FALSE)),file.path(out,"skat-singleton.csv"),row.names=FALSE,quote=FALSE)
# Fitted random-intercept model, independent of GPU eQTL's supplied-covariance GLS.
mixed <- data.frame(group=rep(1:10,each=4),x=rep(c(-1.5,-.5,.5,1.5),10))
mixed$y <- 1.2 + .8*mixed$x + rep(sin(1:10)*1.7,each=4) + cos((1:40)*2.1)*.4
fit <- nlme::lme(y~x,random=~1|group,data=mixed,method="REML",control=nlme::lmeControl(msTol=1e-10))
write.csv(mixed,file.path(out,"mixed-data.csv"),row.names=FALSE,quote=FALSE)
write.csv(data.frame(intercept=nlme::fixef(fit)[1],slope=nlme::fixef(fit)[2],group_variance=as.numeric(nlme::VarCorr(fit)[1,1]),residual_variance=fit$sigma^2,loglik=as.numeric(logLik(fit))),file.path(out,"mixed-reference.csv"),row.names=FALSE,quote=FALSE)
tails <- expand.grid(df=c(1,2,5,30,1000,1e6), t=c(0,1e-8,0.1,1,2,8,30,1e3,1e10))
tails$logp <- pt(tails$t,tails$df,lower.tail=FALSE,log.p=TRUE)
write.csv(tails,file.path(out,"t-tail.csv"),row.names=FALSE,quote=FALSE)
.libPaths(c("target/audit-r-lib",.libPaths()))
library(CompQuadForm)
mix <- data.frame()
for (lambda in list(c(1),c(2,2),c(.5,2,7),c(.1,.5,1,2,4),c(.01,1,10))) {
  for (multiple in c(.1,.5,1,2,5,10)) {
    q <- sum(lambda)*multiple
    if(length(unique(lambda))==1) reference <- pchisq(q/lambda[1],length(lambda),lower.tail=FALSE) else {
      result <- davies(q,lambda,acc=1e-7,lim=1000000)
      stopifnot(result$ifault==0, result$Qq>=0, result$Qq<=1)
      reference <- result$Qq
    }
    mix <- rbind(mix,data.frame(lambda=paste(lambda,collapse=";"),q,p=reference))
  }
}
write.csv(mix,file.path(out,"mixtures.csv"),row.names=FALSE,quote=FALSE)
writeLines(sub("[[:blank:]]+$", "", capture.output(sessionInfo())),file.path(out,"session-info.txt"))
