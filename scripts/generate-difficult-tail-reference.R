# Independent R/CompQuadForm references; all spectra are synthetic.
extra <- Sys.getenv("EQTL_R_LIB")
if (nzchar(extra)) .libPaths(c(extra, .libPaths()))
stopifnot(requireNamespace("CompQuadForm", quietly=TRUE))
stopifnot(packageVersion("CompQuadForm")==package_version("1.4.4"))
spectra <- list(c(1e-12,1,10),c(1e-8,.01,.1,1),
                10^seq(-12,0,length.out=8),10^seq(-12,0,length.out=20),
                10^seq(-12,0,length.out=80),c(rep(1e-12,64),1,2,4),
                c(rep(1e-4,100),rep(1,20)))
rows <- list()
for (lambda in spectra) for (factor in c(.01,.1,1,3,6)) {
  q <- sum(lambda)*factor
  reference <- CompQuadForm::davies(q,lambda,lim=1000000,acc=1e-12)
  if (reference$ifault!=0 || reference$Qq<1e-9) next
  check <- CompQuadForm::imhof(q,lambda,epsabs=1e-11,epsrel=1e-11,limit=100000)
  stopifnot(abs(reference$Qq-check$Qq)<max(3e-10,5*check$abserr))
  rows[[length(rows)+1]] <- data.frame(q=sprintf("%.17g",q),
    eigenvalues=paste(sprintf("%.17g",lambda),collapse=";"),
    survival=sprintf("%.17g",reference$Qq),source="CompQuadForm-1.4.4-Davies-checked-Imhof")
}
out <- "test/resources/difficult-tail-reference.csv"
write.table(do.call(rbind,rows),out,sep=",",row.names=FALSE,quote=FALSE)
cat(length(rows),"independently cross-checked difficult spectra written to",out,"\n")
