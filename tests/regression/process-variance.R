# Run: Rscript --vanilla tests/regression/process-variance.R .
# Requires rjags/JAGS. Tests the actual R generator and its emitted JAGS expression.
args <- commandArgs(trailingOnly=TRUE)
root <- if (length(args)) args[1] else '.'
if (!requireNamespace('rjags',quietly=TRUE)) stop('rjags and JAGS are required; no skip permitted')
prod <- new.env(parent=globalenv())
sys.source(file.path(root,'R','write_JAGS_model.R'),envir=prod)
emit <- function(type='raw',n.iso=2,resid=FALSE) {
  mix <- list(N=if (resid) 2 else 1,n.fe=0,n.re=0,n.iso=n.iso,n.effects=0,
              n.ce=0,FAC=list(),fac_nested=FALSE,factors=NULL,fere=FALSE)
  src <- list(data_type=type,by_factor=NA,conc_dep=FALSE)
  p <- tempfile(fileext='.txt');on.exit(unlink(p))
  prod$write_JAGS_model(p,resid,TRUE,mix,src)
  readLines(p)
}
raw <- emit()
passed <- 0L; failed <- 0L
check <- function(label,expr) {
  tryCatch({force(expr);passed <<- passed+1L;cat('PASS:',label,'\n')},
           error=function(e) {failed <<- failed+1L;cat('FAIL:',label,conditionMessage(e),'\n')})
}
expected <- 'inprod(src_cov[,iso,iso],p2[i,]) + inprod(frac_sig2[,iso],p2[i,])'
check('source covariance parameterization unchanged', {
  stopifnot(any(grepl('src_cov[src,,] <- src_var[src,,] %*% rho[src,,] %*% src_var[src,,]',raw,fixed=TRUE)))
})
check('two-isotope raw process variance uses fitted covariance', {
  stopifnot(any(grepl(expected,raw,fixed=TRUE)))
})
check('one-isotope raw process variance uses fitted covariance', {
  stopifnot(any(grepl(expected,emit(n.iso=1),fixed=TRUE)))
})
check('mean/SD source variance representation unchanged', {
  stopifnot(any(grepl('inprod(1/src_tau[,iso],p2[i,])',emit(type='means'),fixed=TRUE)))
})
check('residual-times-process covariance use unchanged', {
  stopifnot(any(grepl('inprod(src_cov[,i,j],p2[ind,])',emit(resid=TRUE),fixed=TRUE)))
})
check('JAGS evaluates the emitted expression to 50, not 5', {
  line <- raw[grepl('process.var[iso,i] <-',raw,fixed=TRUE)]
  stopifnot(length(line)==1L)
  model <- paste('model {',
    'for (src in 1:2) {',
    'src_var[src,1:2,1:2] <- inverse(src_tau[src,1:2,1:2]);',
    'src_cov[src,1:2,1:2] <- src_var[src,1:2,1:2] %*% rho[src,1:2,1:2] %*% src_var[src,1:2,1:2];',
    '}',
    'for (i in 1:1) { for (iso in 1:2) {',line,'} }',
    'dummy ~ dnorm(0,1)', '}',sep='\n')
  tau <- array(0,c(2,2,2));rho <- tau
  for (s in 1:2) {tau[s,,] <- diag(c(.1,.1));rho[s,,] <- diag(2)}
  con <- textConnection(model)
  j <- tryCatch(rjags::jags.model(con,
    data=list(src_tau=tau,rho=rho,p2=matrix(c(.25,.25),1),frac_sig2=matrix(0,2,2)),
    inits=list(.RNG.name='base::Wichmann-Hill',.RNG.seed=20260906),
    n.chains=1,n.adapt=0,quiet=TRUE),finally=close(con))
  draws <- as.matrix(rjags::coda.samples(j,'process.var',n.iter=2,progress.bar='none'))
  stopifnot(all(abs(draws-50)<1e-10))
})
cat(sprintf('PACK_RESULT={"tests":%d,"failures":%d,"errors":0,"skipped":0}\n',passed+failed,failed))
quit(status=if (failed) 1L else 0L)
