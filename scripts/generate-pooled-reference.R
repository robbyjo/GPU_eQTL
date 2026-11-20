# Copyright 2026 Roby Joehanes; GNU GPL version 3.
# Independent base-R dense covariance reference, not the production spectral algebra.
args <- commandArgs(trailingOnly=TRUE)
out <- if(length(args)) args[1] else 'test/resources/pooled-reference'
dir.create(out, recursive=TRUE, showWarnings=FALSE)
n <- 96L; i <- 0:(n-1L)
group <- rep(c('FHS','WHI','CARDIA'), c(48,24,24))
x <- sin(i*.71)+cos(i*.19)
X <- cbind(Intercept=1,WHI=as.numeric(group=='WHI'),CARDIA=as.numeric(group=='CARDIA'),x=x,
 FHS_batch=as.numeric(group=='FHS' & i%%8>=4),WHI_x2=as.numeric(group=='WHI' & i%%3==0),
 CARDIA_year2=as.numeric(group=='CARDIA' & i%%7==0))
A <- matrix(0,n,n)
for(f in 0:11){j <- (4*f+1):(4*f+4);A[j,j]<-.5;diag(A)[j]<-1}
y <- as.vector(X %*% c(.3,.2,-.1,.07,.08,-.09,.12)) +
 ifelse(group=='FHS', .7*sin(floor(i/4)*1.3)+.32*cos(i*2.1),
 ifelse(group=='WHI',.8*sin(i*1.71),.23*cos(i*1.43)))
g <- (i*7L+i%/%5L)%%3L
data <- data.frame(ID=paste0('s',i),Group=group,y=y,g=g,X[,-1,drop=FALSE],check.names=FALSE)
write.csv(data,file.path(out,'data.csv'),row.names=FALSE,quote=FALSE)
ped <- c('ID\tSIRE\tDAM',unlist(lapply(0:11,function(f)c(paste0('f',f,'\t0\t0'),paste0('m',f,'\t0\t0'),
 sapply(0:3,function(k)paste0('s',f*4+k,'\tf',f,'\tm',f))))))
writeLines(ped,file.path(out,'pedigree.tsv'))
results <- list()
for(case in c('full','missing')) {
 keep <- if(case=='full') seq_len(n) else setdiff(seq_len(n),c(3,8,19,42,53,81,94))
 xx<-X[keep,,drop=FALSE]; yy<-y[keep];gg<-g[keep];nn<-length(keep)
 kernels<-list(A[keep,keep],diag(as.numeric(group[keep]=='CARDIA')),diag(as.numeric(group[keep]=='FHS')),diag(as.numeric(group[keep]=='WHI')))
 scale<-sum(qr.resid(qr(xx),yy)^2)/(nn-ncol(xx)); lo<-log(scale*1e-10);hi<-log(scale*1e10)
 evaluate<-function(theta,gradient=FALSE,return_state=FALSE){
  v<-Reduce('+',Map(function(k,t)k*exp(t),kernels,theta))
  ch<-chol(v); w<-chol2inv(ch);wx<-w%*%xx;m<-crossprod(xx,wx);mi<-chol2inv(chol(m))
  p<-w-wx%*%mi%*%t(wx);py<-as.vector(p%*%yy)
  value<-.5*((nn-ncol(xx))*log(2*pi)+2*sum(log(diag(ch)))+as.numeric(determinant(m,logarithm=TRUE)$modulus)+sum(yy*py))
  grad<-sapply(seq_along(kernels),function(j).5*exp(theta[j])*(sum(p*t(kernels[[j]]))-sum(py*(kernels[[j]]%*%py))))
  if(return_state)return(list(p=p,py=py,value=value))
  if(gradient)grad else value
 }
 fits<-lapply(c(.02,.2,.7),function(a)nlminb(log(scale*c(a,.08,.4,.8)),evaluate,function(t)evaluate(t,TRUE),lower=rep(lo,4),upper=rep(hi,4),control=list(eval.max=2000,iter.max=1000,rel.tol=1e-13,x.tol=1e-11)))
 fit<-fits[[which.min(sapply(fits,`[[`,'objective'))]];state<-evaluate(fit$par,return_state=TRUE)
 U<-sum(gg*state$py);info<-as.numeric(crossprod(gg,state$p%*%gg));stat<-U^2/info
 results[[case]]<-data.frame(Case=case,N=nn,Genetic=exp(fit$par[1]),CARDIA=exp(fit$par[2]),FHS=exp(fit$par[3]),WHI=exp(fit$par[4]),
 Effect=U/info,SE=1/sqrt(info),Statistic=stat,log10P=pchisq(stat,1,lower.tail=FALSE,log.p=TRUE)/log(10),REMLLogLik=-state$value)
}
write.table(do.call(rbind,results),file.path(out,'expected.csv'),sep=',',row.names=FALSE,col.names=TRUE,quote=FALSE)
writeLines(sub('[[:space:]]+$','',capture.output(sessionInfo())),file.path(out,'session-info.txt'))
